#include <Rcpp.h>
#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "../Cpp/Routines/C++/RandomNumbers/random.h"

using namespace Rcpp;
using namespace std;

#define GET_VARIABLE_NAME(Variable) (#Variable)

enum client { resident, visitor, absence };																		
// clients can be resident, visitors, or be absent
enum learPar { alphaPar, gammaPar, netaPar , alphathPar};

enum learnScenario {nature, experiment, marketExperiment, ExtendedMarket};


class agent													// Learning agent
{
public:
  agent(double alphaI, double gammaI, double negRewardI, 
        double alphathI, double initVal);
  // constructor providing values for the learning parameters
  ~agent();																
  // destructor not really necessary
  void update();
  // function that updates the value of state-action pairs according to 
  //current reward and estimates of future values
  void act(client newOptions[], int &idNewOptions, double &VisProbLeav, 
           double ResProbLeav, double VisReward, double ResReward, double inbr, 
           double outbr,  learnScenario scenario);
  // function where the agent takes the action, gets reward, see new 
  //state and chooses future action
  double getLearnPar(learPar parameter);
  // function to access the learning parameters
  void checkChoice();
  // Check that the choice taken is among one of the options, 
  //otherwise trigger an error
  void rebirth(double initVal);		
  // Function to reset private variables in an individual
  bool getChoice(bool time) {
    if (!time) return choiceT;
    else return choiceT1;
  }
  void ObtainReward(double &ResReward, double &VisReward);
  // Get reward
  void getNewOptions(client newOptions[], int &idNewOptions, 
                     double &VisProbLeav, double &ResProbLeav,  
                     double &inbr, double &outbr, learnScenario &scenario);
  // Function to get new clients in the station, when in a natural environment
  void getExternalOptions(client newOptions[], int &idNewOptions, 
                          double &inbr, double &outbr);		
  // After unattended clients leave or stay, get new clients
  void getExperimentalOptions();
  // Get new clients in the experimental setting
  void getMarketExperiment();
  // Get new clients in the experimental setting of Olle's model
  void getExtenededMarket();
  // Get new clients in the experimental setting of Noa's experiment
  double logist();
  int mapOptionsDP(client options[], int &choice);			
  bool isRVoption(int time);
  // Function to know if there is an RV option in the 
  // current t=0 or future t=1 option
  client cleanOptionsT[2];	// current cleaning options time = t
  client cleanOptionsT1[2];	// future cleaning options  time = t+1
  void choice();
  // Function to make a choice
  virtual int mapOptions(client options[], int &choice)=0;
  // function that maps state action pairs to indexes in the array 'values' 
  //where values are stored
  virtual void updateThet(int curState) = 0;
  // function to update the policy parameter Theta
  int numEst;
  // Number of estimates characterizing behavioural options 9 for FAA
  int countExp;
protected:
  double values[6];																							
  // array storing the estimated values of states 
  double delta;
  double piV;
  double theta[2]; // policy parameters
  int DPid;
  int choiceT;// current choice 
  int choiceT1;// future choice
  double alpha;// speed of learning for estimated values
  double alphath; // speed of learning for policy parameter
  double gamma;// importance of future rewards
  bool neta=1;
  double negRew_const;
  // Weight of the negative reward in the total reward obtained by an agent
  double currentReward; // reward given by current state action pair
  double cumulReward;	// Cumulative reward
  int age;
  double negRew_curr; // size of the negative reward in the current time step
  int currState, nextState;
};

void wait_for_returnRcpp(string strTOp="return"){
  Environment base = Environment("package:base");
  Function readline = base["readline"];
  Function as_numeric = base["as.numeric"];
  std::cout << strTOp;
  int drink = as<int>(as_numeric(readline("> ")));
}

// Members of agent class

agent::agent(double alphaI = 0.01, double gammaI = 0.5, 
             double negRewardI = 0, double alphathI = 0.01,
             double initVal = 0){
  // parameterized constructor with defaults
  theta[0] = 0, theta[1] = 0;
  numEst = 6;
  delta = 0;
  for (int i = 0; i < numEst; i++) { values[i] = 1+initVal; }
  // assigned educated initial values. Reward + guess of 
  //the expected future reward
  values[5] -= 1;
  // Value of absence starts with reward of 0
  piV = logist();
  alpha = alphaI, gamma = gammaI, alphath = alphathI;
  negRew_const = negRewardI;
  cleanOptionsT[0] = absence, cleanOptionsT[1] = absence, choiceT = 0;
  cleanOptionsT1[0] = absence, cleanOptionsT1[1] = absence, choiceT1 = 0;
  currentReward = 0, cumulReward = 0; negRew_curr = 0;
  age = 0;
  countExp = 2;
}

void agent::rebirth(double initVal = 0){
  currentReward = 0;
  cumulReward = 0;
  choiceT = 0, choiceT1 = 0;
  cleanOptionsT[0] = absence, cleanOptionsT[1] = absence;
  cleanOptionsT1[0] = absence, cleanOptionsT1[1] = absence;
  age = 0;
  for (int i = 0; i < numEst; i++) { values[i] = 1 + initVal; }
  values[numEst-1] -= 1;
  piV = logist();
  delta = 0;
  theta[0] = 0, theta[1] = 0;
  countExp = 2;
}

agent::~agent() {}		// Destructor

void agent::checkChoice(){
  if (choiceT > 1 )  {
    error("agent::act", "choice is not among the options");
  }
}

double agent::getLearnPar(learPar parameter)
{
  switch (parameter)
  {
  case alphaPar:return(alpha);
    break;
  case gammaPar:
    return(gamma);
    break;
  case netaPar:return(neta);
    break;
  case alphathPar:return(alphath);
    break;
  default:error("agent:getlearnPar",
                "asking for a parameter that does not exist");
  return 0;
  break;
  }
}

void agent::ObtainReward(double &ResReward, double &VisReward){
  if (cleanOptionsT[choiceT] == resident) { 
    currentReward = ResReward, cumulReward += ResReward; 
  }					// Obtain reward if the choice is a resident
  else if (cleanOptionsT[choiceT] == visitor) { 
    currentReward = VisReward, cumulReward += VisReward; 
  }					// Obtain reward if the choice is a visitor
  else { currentReward = 0, cumulReward += 0; }
  // No reward if there is no client in the choice
}

// Function to get new clients in the station, when in a natural environment
void agent::getNewOptions(client newOptions[], int &idNewOptions, 
                          double &VisProbLeav, double &ResProbLeav, double &inbr, double &outbr, 
                          learnScenario &scenario){
  if (choiceT == 0)		// Define the behaviour of the unattended client
  {
    if (cleanOptionsT[1] == resident)
    {
      if (rnd::uniform() > ResProbLeav) { 
        cleanOptionsT1[0] = cleanOptionsT[1], negRew_curr = 0; 
      }								
      // if the unttended client is a resident, it leaves with probability ResPropLeave
      else { negRew_curr = negRew_const; }
    }
    else if (cleanOptionsT[1] == visitor)
    {
      if (rnd::uniform() > VisProbLeav) { 
        cleanOptionsT1[0] = cleanOptionsT[1], negRew_curr = 0; 
      }
      // if the unttended client is a visitor, it leaves with probability VisPropLeave
      else { 
        negRew_curr = negRew_const; 
      }
    }
    else { negRew_curr = 0; }
  }
  else
  {
    if (cleanOptionsT[0] == resident)
    {
      if (rnd::uniform() > ResProbLeav) { 
        cleanOptionsT1[0] = cleanOptionsT[0], negRew_curr = 0;
      }	
      // if the unattended client is a resident, it leaves with probability ResPropLeave
      else { negRew_curr = negRew_const; }
    }
    else if (cleanOptionsT[0] == visitor)
    {
      if (rnd::uniform() > VisProbLeav) { 
        cleanOptionsT1[0] = cleanOptionsT[0], negRew_curr = 0; }
      // if the unattended client is a visitor, it leaves with probability VisPropLeave
      else { negRew_curr = negRew_const; }
    }
    else { negRew_curr = 0; }
  }
  switch (scenario) {
  case nature: 
    getExternalOptions(newOptions, idNewOptions, inbr, outbr);
    break;
  case experiment: 
    getExperimentalOptions();
    break;
  case marketExperiment: 
    getMarketExperiment();
    break;
  case ExtendedMarket: 
    getExtenededMarket();
    break;
  default:cout << "unkown scenario!!" << endl;
  error("agent:getNewOptions",
        "unkown scenario");
  break;
  }
}

void agent::getExternalOptions(client newOptions[], int &idNewOptions, 
                               double &inbr, double &outbr)
{
  if (cleanOptionsT1[0] == absence){
    // If none of the clients stayed from the previous interaction
    cleanOptionsT1[0] = newOptions[idNewOptions], ++idNewOptions;
    if (cleanOptionsT1[0] == absence){
      // If the first draw does not yield a client
      cleanOptionsT1[1] = newOptions[idNewOptions], ++idNewOptions;
      return;
    }
  }
  if (cleanOptionsT1[0] != absence){
    // Fill the second option depending on the first option
    double probs[3] = { (1 - inbr)*(1 - outbr) + inbr*outbr , 0, 0 };	
    // Define probabilities depending on parameters
    probs[1] = probs[0] + inbr*(1 - outbr);
    // First prob is of a random option	
    probs[2] = probs[1] + outbr*(1 - inbr);	
    // Second and third homophily, and heterophily respectively
    if (probs[2] != 1) error("agent:getExternalOptions", 
        "probability does not sum up to 1");
    double rand = rnd::uniform();
    if (probs[0] > rand) {	
      cleanOptionsT1[1] = newOptions[idNewOptions], ++idNewOptions;
    }						// Random
    else if (probs[1] > rand)	cleanOptionsT1[1] = cleanOptionsT1[0];												
    // homophily
    else					// heterophily
    {
      if (cleanOptionsT1[0] == resident) { cleanOptionsT1[1] = visitor; }
      else { cleanOptionsT1[1] = resident; }
    }
  }
}

void agent::getExperimentalOptions() {
  // Get new options in an experimental setting
  if (cleanOptionsT[0] == resident && cleanOptionsT[1] == visitor) {
    return;	
  }	// Every other option is a Resident-Visitor
  else {
    cleanOptionsT1[0] = resident;
    cleanOptionsT1[1] = visitor;
    return;
  }
}

void agent::getMarketExperiment() {
  // Get new options in an experimental setting of Olle's models
  if (countExp==0){
    //(cleanOptionsT[0] == resident && cleanOptionsT[1] == visitor) {
    ++countExp;
    return;
  }	// Every other option is a Resident-Visitor
  else if (countExp==1){
    /*((cleanOptionsT[0] == resident || cleanOptionsT[0] == visitor)
     && cleanOptionsT[1] == absence) {*/
    cleanOptionsT1[0] = absence;
    cleanOptionsT1[1] = absence;
    ++countExp;
    return;
  }
  else {
    cleanOptionsT1[0] = resident;
    cleanOptionsT1[1] = visitor;
    countExp = 0;
  }
}

void agent::getExtenededMarket() {
  // Get new options in the experimental setting of Noa's experiment
  if (countExp==0){
    //(cleanOptionsT[0] == resident && cleanOptionsT[1] == visitor) {
    ++countExp;
    return;
  }	// Every other option is a Resident-Visitor
  else if (countExp==1){
    /*((cleanOptionsT[0] == resident || cleanOptionsT[0] == visitor)
     && cleanOptionsT[1] == absence) {*/
    cleanOptionsT1[0] = absence;
    cleanOptionsT1[1] = absence;
    ++countExp;
    return;
  }
  else {
    countExp = 0;
    double rand = rnd::uniform();
    if (rand < 0.5) {
      cleanOptionsT1[0] = resident;
      cleanOptionsT1[1] = visitor;
    }
    else if (rand<0.75) {
      cleanOptionsT1[0] = resident;
      cleanOptionsT1[1] = resident;
    }
    else {
      cleanOptionsT1[0] = visitor;
      cleanOptionsT1[1] = visitor;
    }
  }
}

void agent::act(client newOptions[], int &idNewOptions, double &VisProbLeav, 
                double ResProbLeav, double VisReward, double ResReward, double inbr,
                double outbr, learnScenario scenario){
  // taking action, obatining reward, seeing new state, choosing future action
  ++age;
  // new time step
  cleanOptionsT[0] = cleanOptionsT1[0], cleanOptionsT[1] = cleanOptionsT1[1];
  // Future state becomes current state
  choiceT = choiceT1;
  // Future action becomes current action
  checkChoice();
  // Check that the choice is among the options
  cleanOptionsT1[0] = absence, cleanOptionsT1[1] = absence;
  // Future state is unknown
  choiceT1 = 2;
  ObtainReward(VisReward,ResReward);
  getNewOptions(newOptions, idNewOptions, VisProbLeav, ResProbLeav,
                inbr, outbr, scenario);
  choice();
}

void agent::update(){																								
  // change estimated value according to TD error
  // change policy parameter according to TD error
  // why did this following lines changed while developing the 
  currState = mapOptions(cleanOptionsT,choiceT);
  nextState = mapOptions(cleanOptionsT1, choiceT);
  delta = currentReward -
    negRew_curr*neta + gamma*values[nextState] - values[currState];
  // construct the TD error
  values[currState] += alpha*delta;
  // update value
  updateThet(currState);
}

double agent::logist() { return (1 / (1 + exp(-(theta[0]-theta[1]))));}

int agent::mapOptionsDP(client options[], int &choice){
  int state;
  if (options[0] == absence || options[1] == absence)	{
    // One of the options is empty
    if (options[0] == resident || options[1] == resident){					
      // the other one is a resident
      state = 2;                                                    // R0
    }
    else if (options[0] == visitor || options[1] == visitor){
      // the other one is a visitor
      state = 1;                                                    // V0
    }
    else { state = 5; }				                                  // 00
  }
  else if (options[0] == resident || options[1] == resident){
    // Both options have clients and one of them is a resident
    if (options[0] == visitor || options[1] == visitor){
      // the other one is a visitor
      state = 0;                                                    // RV
    }
    else { state = 3; }		                                          // RR
  }
  else { state = 4; }			                                          // VV
  return state;
}

void agent::choice() {
  if (cleanOptionsT1[0] == absence || cleanOptionsT1[1] == absence) {
    // if there is an absence choose the client
    bool presence = cleanOptionsT1[0] == absence;
    choiceT1 = presence;
  }
  else if (cleanOptionsT1[0] != cleanOptionsT1[1]) {
    // if clients are different use policy (logist)
    bool visit = rnd::bernoulli(piV);
    if (cleanOptionsT1[1] == visitor) {
      choiceT1 = visit;
    } else {
      choiceT1 = !visit;
    }
  }
  else {
    // if the clients are the same - choose randomly
    choiceT1 = rnd::bernoulli();
  }
}

bool agent::isRVoption(int time) {
  if (time == 0) {
    return(cleanOptionsT[0] + cleanOptionsT[1] == 1);
  }
  else {
    return(cleanOptionsT1[0] + cleanOptionsT1[1] == 1);
  }
}

class FAATyp1 :public agent{			// Fully Aware Agent (FAA)			
public:
  FAATyp1(double alphaI, double gammaI, double NegRewI, 
          double alphaThI, double initVal=1):agent(alphaI, gammaI, NegRewI, 
          alphaThI, initVal){
  }
  virtual int mapOptions(client options[], int &choice){
    return(mapOptionsDP(options, choice));
  }

  virtual void updateThet(int curState) {
    if (curState == 0) {
      if (cleanOptionsT[choiceT] == visitor) {
        theta[0] += alphath*delta*(1 - piV);
        theta[1] -= alphath*delta*(1 - piV);
      } else {
        theta[0] -= alphath*delta*piV;
        theta[1] += alphath*delta*piV;
      }
      piV = logist();
    }
  }
};

class PAATyp1 :public agent{				// Partially Aware Agent (PAA)	
public:
  PAATyp1(double alphaI, double gammaI, double netaI, 
          double alphaThI, double initVal, double alphaThNchI):agent(alphaI, gammaI,  
          netaI, alphaThI,initVal){
    alphaThNch = alphaThNchI;
    numEst = 3;
    values[2] -= 1;
  }
  
  int mapOptions(client options[], int &choice){
    if (options[choice] == resident) { return (0); }
    else if (options[choice] == visitor) { return(1); }
    else { return(2); }
    return(options[choice]); 
  }
  virtual void updateThet(int curStatAct) {
    if (curStatAct < 2) {
      //int notchoice = (choiceT == 0);
      if (curStatAct == 1) {
        if (cleanOptionsT[0] == cleanOptionsT[1]) {
          theta[0] += alphaThNch*alphath*delta*(1 - piV);
        }
        else {
          theta[0] += alphath*delta*(1 - piV);
        }
        
      }
      else {
        if (cleanOptionsT[0] == cleanOptionsT[1]) {
          theta[1] += alphaThNch*alphath*delta*piV;
        }
        else {
          theta[1] += alphath*delta*piV;
        }
        
      }
      piV = logist();
    }
  }
private:
  double alphaThNch;
};


Rcpp::DataFrame abs2rel_abund(Rcpp::NumericVector abund_clean,
                              Rcpp::NumericVector abund_resid,
                              Rcpp::NumericVector abund_visitors,
                              Rcpp::List param) {
  double totAbundance_inv, inv_tot_abund_client;
  Rcpp::NumericVector rel_abund_clean;
  Rcpp::NumericVector rel_abund_resid;
  Rcpp::NumericVector rel_abund_visitors;
  for (int d_point = 0; d_point < abund_clean.size(); ++d_point){
    totAbundance_inv = 1/(double(param["scaleConst"])*abund_clean[d_point] + 
      abund_resid[d_point] +
      abund_visitors[d_point]);
    inv_tot_abund_client = 1 / (abund_resid[d_point] +
      abund_visitors[d_point]);
    rel_abund_clean.push_back(
      double(param["scaleConst"])*abund_clean[d_point]*totAbundance_inv);
    rel_abund_resid.push_back((1- rel_abund_clean[d_point])*
      abund_resid[d_point]*inv_tot_abund_client);
    rel_abund_visitors.push_back((1 - rel_abund_clean[d_point])*
      abund_visitors[d_point]*inv_tot_abund_client);
  }
  return Rcpp::DataFrame::create(_("rel_abund_clean")=rel_abund_clean,
                                 _("rel_abund_resid")=rel_abund_resid,
                                 _("rel_abund_visitors")=rel_abund_visitors);
}

void checkGroups(Rcpp::IntegerVector &group,bool groups) {
  if (!groups) {
    for (int i = 0; i < group.size(); ++i)
      group[i] = 0;
  }
  return;
}

bool setAgent(Rcpp::List currPars, Rcpp::List sim_param, int group,
              double rel_abund_clean,double prob_Vis_Leav) {
  
  double prob_F;
  Rcpp::IntegerVector probFAAs = currPars["probFAA"];
  switch (int(sim_param["agentScen"])){
  case 0:
    prob_F = probFAAs[group];
    break;
  case 1:
    prob_F = 1 / (1 + exp(-(double(currPars["interpReg"])
                              + double(currPars["slopRegRelAC"])*rel_abund_clean
                              + double(currPars["slopRegPVL"])*prob_Vis_Leav)));
  default:
    prob_F =probFAAs[group];
  break;
  }
  return rnd::bernoulli(prob_F);
}

void draw(client trainingSet[], int rounds, double probRes, double probVis){					
  // In a natural setting draw clients according to their abundance
  double cumProbs[3] = { probRes, probRes + probVis, 1 };
  double rndNum;
  for (int i = 0; i < rounds * 2; i++){
    rndNum = rnd::uniform();
    if (rndNum < cumProbs[0]) { trainingSet[i] = resident; }
    else if (rndNum < cumProbs[1]) { trainingSet[i] = visitor; }
    else { trainingSet[i] = absence; }
  }
}

// [[Rcpp::export]]
Rcpp::DataFrame do_simulation(
    Rcpp::DataFrame emp_data, Rcpp::List focal_param, Rcpp::List sim_param) {
  client clientSet[int(sim_param["totRounds"]) * 2];
  // array to store the queue of clients
  int idClientSet;
  // counter for the client queue
  Rcpp::NumericVector gammas = focal_param["gamma"];
  Rcpp::NumericVector negRewards = focal_param["negReward"];
  agent* cleaners[2][2];
  PAATyp1 PAAcleaner_NG(focal_param["alphaC"], gammas[0],
                        negRewards[0], focal_param["alphaA"],
                                                  1.0, 0.0);
  cleaners[0][0] = &PAAcleaner_NG;
  FAATyp1 FAAcleaner_NG (focal_param["alphaC"], gammas[0],
                         negRewards[0], focal_param["alphaA"]);
  cleaners[1][0] = &FAAcleaner_NG;
  if (sim_param["Group"]) {
    PAATyp1 PAAcleaner_G (focal_param["alphaC"], gammas[1],
                          negRewards[1], focal_param["alphaA"], 1.0, 0.0);
    cleaners[0][1] = &PAAcleaner_G;
    FAATyp1 FAAcleaner_G (focal_param["alphaC"], gammas[1],
                          negRewards[2], focal_param["alphaA"]);
    cleaners[1][1] = &FAAcleaner_G;
  }
  double VisPref, init;
  int countRVopt;
  Rcpp::DataFrame relAbunds = abs2rel_abund(emp_data["abund_clean"],
                                            emp_data["abund_resid"],
                                                    emp_data["abund_visitors"],
                                                            focal_param);
  Rcpp::IntegerVector group = emp_data["group"];
  Rcpp::NumericVector rel_abund_clean = relAbunds["rel_abund_clean"];
  Rcpp::NumericVector prob_Vis_Leav = emp_data["prob_Vis_Leav"];
  Rcpp::NumericVector rel_abund_resid = relAbunds["rel_abund_resid"];
  Rcpp::NumericVector rel_abund_visitors = relAbunds["rel_abund_visitors"];
  Rcpp::LogicalVector agent;
  Rcpp::NumericVector marketPred;
  checkGroups(group,sim_param["Group"]);
  for (int id_data_point = 0; id_data_point < 5; ++id_data_point) {
    // emp_data.nrows(); ++id_data_point) {
    agent.push_back(setAgent(focal_param,sim_param,
                             group[id_data_point],
                                  rel_abund_clean[id_data_point],
                                                 prob_Vis_Leav[id_data_point]));
    init =   gammas[group[id_data_point]]*
      (1 - pow(1 -
      rel_abund_resid[id_data_point] -
      rel_abund_visitors[id_data_point], 2)) /  (1 - gammas[group[id_data_point]]);
    draw(clientSet, int(sim_param["totRounds"]),
         rel_abund_resid[id_data_point],
                        rel_abund_visitors[id_data_point]);
    cleaners[agent[id_data_point]][group[id_data_point]]->rebirth(init);
    VisPref = 0, countRVopt = 0,  idClientSet = 0;
    // Loop through the learning rounds
    for (int trial = 0; trial < int(sim_param["totRounds"]); ++trial) {
      cleaners[agent[id_data_point]][group[id_data_point]]->
        act(clientSet, idClientSet, prob_Vis_Leav[id_data_point],
            double(sim_param["ResProbLeav"]), double(sim_param["VisReward"]),
            sim_param["ResReward"], sim_param["inbr"], sim_param["outbr"],
            learnScenario(int(sim_param["scenario"])));
      cleaners[agent[id_data_point]][group[id_data_point]]->update();
      if (trial > int(sim_param["totRounds"]) * float(sim_param["propfullPrint"])) {
        if (cleaners[agent[id_data_point]][group[id_data_point]]->isRVoption(0)) {
          ++countRVopt;
          if (cleaners[agent[id_data_point]][group[id_data_point]]->
          cleanOptionsT[cleaners[agent[id_data_point]]
                [group[id_data_point]]->getChoice(0)] == visitor)
            ++VisPref;
          }
        }
      }
      if (countRVopt == 0) marketPred.push_back(0.5);
      else marketPred.push_back( VisPref / countRVopt);
  }
  

  
  // cleaners[0][0]->rebirth(init);
  // VisPref = 0, countRVopt = 0;
  // cleaners[0][0]->
  //   act(clientSet, idClientSet, VisProbLeavDebug,
  //   double(sim_param["ResProbLeav"]), double(sim_param["VisReward"]),
  //   sim_param["ResReward"], sim_param["inbr"], sim_param["outbr"],
  //   learnScenario(int(sim_param["scenario"])));
  // cleaners[0][0]->update();
  // if (trial > int(sim_param["totRounds"]) * float(sim_param["propfullPrint"])) {
  //           if (cleaners[0][0]->
  //           cleanOptionsT[cleaners[0][0]->getChoice(0)] == visitor)
  //             ++VisPref;
  // }
  return(relAbunds);
  
  // Loop through the learning rounds
  // for (int trial = 0; trial < int(sim_param["totRounds"]); ++trial) {
  //     
  //     }
  //     if (countRVopt == 0) marketPred.push_back(0.5);
  //     else marketPred.push_back( VisPref / countRVopt);
  //     cleaners[agent[id_data_point]][group[id_data_point]]->rebirth();
  //   }
  // }
  // delete[] clientSet;
  // emp_data.push_back(rel_abund_clean,"rel_abund_clean");
  // emp_data.push_back(rel_abund_resid,"rel_abund_resid");
  // emp_data.push_back(rel_abund_visitors,"rel_abund_visitors");
  // emp_data.push_back(group,"group");
  // emp_data.push_back(agent,"agent");
  // emp_data.push_back(marketPred,"marketPred");
  // return emp_data;
}
