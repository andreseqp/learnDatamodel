 // ActCrit.cpp : Defines the entry point for the console application.
//
/*=============================================================================
ActCrit.cpp
===============================================================================

Main file of a project to use data to fit parameter values in a learning model.
The fit is done using a Markov Chain Monte Carlo. Address of the data file
and parameter values for the MCMC arer provided to the executable file through 
a JSON file. The output of the executable is either the MCMC data or a data file 
including predictions from the computational model. 



Written by:

Andr�s E. Qui�ones
Posdoctoral researcher
Behavioural Ecology Group
Institute of Biology
University of Neuch�tel
Switzerland

Start date:
5 April 2017

=============================================================================*/

#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include <iostream>
#include <fstream>
#include <omp.h>
#include "random.h"
#include "json.hpp"       
// Header for reading and using JSON files see https://github.com/nlohmann/json


#define GET_VARIABLE_NAME(Variable) (#Variable)

using namespace std;
using json = nlohmann::json;

// General parameters

// Classes

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
	void printIndData(ofstream &learnSeries, int &seed, double &outbr, 
		double pV, double pR);
	// prints individual data from the learning process
	double getLearnPar(learPar parameter);
	// function to access the learning parameters
	void checkChoice();
	// Check that the choice taken is among one of the options, 
	//otherwise trigger an error
	void rebirth(double initVal);		
	int getstate(bool time) {
		if (!time) return currState;
		else return nextState;
	}
	bool getChoice(bool time) {
		if (!time) return choiceT;
		else return choiceT1;
	}
	// Function to reset private variables in an individual
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
	void ObtainReward(double &ResReward, double &VisReward);
	// Get reward
	double logist();
	int mapOptionsDP(client options[], int &choice);			
	// default function that maps state pairs to indexes in the array 
	//'values' where values are stored works for DPupdate and for 
	//state-action pair NOT for action estimation
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

void agent::rebirth(double initVal = 0)
{
	age = 0;
	cleanOptionsT[0] = absence, cleanOptionsT[1] = absence;
	cleanOptionsT1[0] = absence, cleanOptionsT1[1] = absence;
	choiceT = 0, choiceT1 = 0;
	currentReward = 0;
	cumulReward = 0;
	for (int i = 0; i < numEst; i++) { values[i] = 1 + initVal; }
	values[5] -= 1;
	piV = logist();
	delta = 0;
	theta[0] = 0, theta[1] = 0;
	countExp = 2;
}

agent::~agent() {}		// Destructor

void agent::checkChoice()
{
	if (choiceT > 1 )
	{
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

void agent::ObtainReward(double &ResReward, double &VisReward)
{
	if (cleanOptionsT[choiceT] == resident) { 
		currentReward = ResReward, cumulReward += ResReward; 
	}					// Obtain reward if the choice is a resident
	else if (cleanOptionsT[choiceT] == visitor) { 
		currentReward = VisReward, cumulReward += VisReward; 
	}					// Obtain reward if the choice is a visitor
	else { currentReward = 0, cumulReward += 0; }
	// No reward if there is no client in the choice
}

void agent::getNewOptions(client newOptions[], int &idNewOptions, 
	double &VisProbLeav, double &ResProbLeav, double &inbr, double &outbr, 
	learnScenario &scenario)
{
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
	case nature: getExternalOptions(newOptions, idNewOptions, inbr, outbr);
		break;
	case experiment: getExperimentalOptions();
		break;
	case marketExperiment: getMarketExperiment();
		break;
	case ExtendedMarket: getExtenededMarket();
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
	currState = mapOptions(cleanOptionsT,choiceT);
	nextState = mapOptions(cleanOptionsT1, choiceT);
	delta = currentReward -
		negRew_curr*neta + gamma*values[nextState] - values[currState];
	// construct the TD error
	values[currState] += alpha*delta;
	// update value
	updateThet(currState);
}

void agent::printIndData(ofstream &learnSeries, int &seed, double &outbr, 
	double pV, double pR)
{
	learnSeries << seed << '\t' << age << '\t';
	//cout << seed << '\t' << age << '\t';
	learnSeries << alpha << '\t' << gamma << '\t';
	learnSeries << neta << '\t' << alphath << '\t' << pV << '\t';
	learnSeries << pR << '\t' << theta[0] << '\t';
	learnSeries << theta[1] << '\t' << outbr << '\t';
	learnSeries << cleanOptionsT[0] << '\t' << cleanOptionsT[1] << '\t';
	learnSeries << cleanOptionsT[choiceT] << '\t';
	//cout << cleanOptionsT[0] << '\t' << cleanOptionsT[1] << '\t' << choiceT << '\t';
	learnSeries << currentReward << '\t' << cumulReward << '\t' << negRew_curr << '\t';
	//cout << currentReward << '\t' << cumulReward << '\t';
	for (int j = 0; j < numEst; j++)
	{
		learnSeries << values[j] << '\t';
		//cout << values[j] << '\t';
	}
	learnSeries << endl;
	//cout << endl;
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
		values[3] -= 1;
	}
	
	void rebirth(int initVal=1) {
		rebirth();
		values[3] -= 1;
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

// Functions external to the agent

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

string create_filename(std::string filename,nlohmann::json param){
	// name the file with the parameter specifications
	filename.append("MCMCchain_CL");
	filename.append(itos(param["chain_length"]));
	filename.append("_seed");
	filename.append(itos(param["seed"]));
filename.append(".txt");
return(filename);
}



struct locat_point {
	double rel_abund_clean, rel_abund_visitors, rel_abund_resid, prob_Vis_Leav;
	string site_year;
	int countMarket, totalMarket;
	float market_exp_success;
	double marketPred;
};

struct cleaner_point {
	double abund_clean, abund_visitors, abund_resid, prob_Vis_Leav;
	double rel_abund_clean, rel_abund_visitors, rel_abund_resid;
	string site_year;
	string cleanerID;
	int countVisitor, totalMarket;
	float market_exp_success;
	double marketPred;
	int group=0;
};



struct model_param {
	//model_param(model_param const &obj);
	double alphaC, alphaA, scaleConst;
	double gamma[2], negReward[2];
	model_param &operator= (model_param const &rhs) {
		alphaA = rhs.alphaA;
		alphaC = rhs.alphaC;
		scaleConst = rhs.scaleConst;
		gamma[0] = rhs.gamma[0];
		negReward[0] = rhs.negReward[0];
		gamma[1] = rhs.gamma[1];
		negReward[1] = rhs.negReward[1];
		return *this;
	}
};

//model_param::model_param(model_param const &obj) {
//	alphaA = obj.alphaA;
//	alphaC = obj.alphaC;
//	gamma = obj.gamma;
//	negReward = obj.negReward;
//}

vector <locat_point> read_locData(ifstream& marketData) {
	// open file
	if (!marketData.is_open()) {
		std::cerr << "error: unable to open data file\n";
		wait_for_return();
		exit(EXIT_FAILURE);
	}
	string header;
	getline(marketData, header); // skip header
	vector<locat_point> data_set;
	locat_point input;
	for (;;) { // read data
		// if end of file
		marketData >> input.site_year;
		marketData >> input.rel_abund_clean;
		marketData >> input.rel_abund_visitors;
		marketData >> input.rel_abund_resid;
		marketData >> input.prob_Vis_Leav;
		marketData >> input.market_exp_success;
		marketData >> input.countMarket;
		marketData >> input.totalMarket;
		if (marketData.eof()) 	break;
		data_set.emplace_back(input);
	}
	return(data_set);
}

vector <cleaner_point> read_CleanData(ifstream& marketData) {
	// open file
	if (!marketData.is_open()) {
		std::cerr << "error: unable to open data file\n";
		wait_for_return();
		exit(EXIT_FAILURE);
	}
	string header;
	getline(marketData, header); // skip header
	vector<cleaner_point> data_set;
	cleaner_point input;
	for (;;) { // read data
		marketData >> input.site_year;
		marketData >> input.cleanerID;
		marketData >> input.countVisitor;
		marketData >> input.abund_clean;
		marketData >> input.abund_visitors;
		marketData >> input.abund_resid;
		marketData >> input.prob_Vis_Leav;
		marketData >> input.group;
		if (marketData.eof()) 	break;
		input.totalMarket = 20;
		input.market_exp_success = double(input.countVisitor) /
			double(input.totalMarket);
		data_set.emplace_back(input);
		// if end of file
	}
	return(data_set);
}

void abs2rel_abund(vector<cleaner_point> & emp_data, model_param param) {
	double totAbundance_inv, inv_tot_abund_client;
	for (int d_point = 0; d_point < emp_data.size(); ++d_point){
		totAbundance_inv = 1/(param.scaleConst*emp_data[d_point].abund_clean + emp_data[d_point].abund_resid +
			emp_data[d_point].abund_visitors);
		inv_tot_abund_client = 1 / (emp_data[d_point].abund_resid +
			emp_data[d_point].abund_visitors);
		emp_data[d_point].rel_abund_clean = 
			param.scaleConst*emp_data[d_point].abund_clean*totAbundance_inv;
		emp_data[d_point].rel_abund_resid =
			(1- emp_data[d_point].rel_abund_clean)*emp_data[d_point].abund_resid*inv_tot_abund_client;
		emp_data[d_point].rel_abund_visitors=
			(1 - emp_data[d_point].rel_abund_clean)*emp_data[d_point].abund_visitors*inv_tot_abund_client;
	}
}

//vector<data_point> read_Data(ifstream &marketData) {
//	// open file
//	if (!marketData.is_open()) {
//		std::cerr << "error: unable to open data file\n";
//		wait_for_return();
//		exit(EXIT_FAILURE);
//	}
//	string header;
//	getline(marketData,header); // skip header
//	vector<data_point> data_set; 
//	data_point input;
//	for (;;) { // read data
//		marketData >> input.rel_abund_clean;
//		marketData >> input.rel_abund_visitors;
//		marketData >> input.rel_abund_resid;
//		marketData >> input.prob_Vis_Leav;
//		marketData >> input.market_exp_success;
//		data_set.emplace_back(input);
//		// if end of file
//		if (marketData.eof()) 	break;
//	}
//	return(data_set);
//}

void initializeChainFile(ofstream &chainOutput,nlohmann::json param){
	std::string namedir = param["folder"];
	string IndFile = create_filename(namedir, param);
	chainOutput.open(IndFile.c_str());
	chainOutput << "iteration	" << "alpha_actor	" << "alpha_critic	" <<
		 "gamma	" << "negReward	";
	if (bool(param["Group"]))
		chainOutput	<< "gamma.1	" << "negReward.1	";
	chainOutput  << "scaleConst	" << "fit	" << "ratio" << endl;
}

void initializeRoundFile(ofstream& RoundOutput, nlohmann::json param,
	model_param focals) {
	// File to record the predictions and data for one parameter combination
	std::string namedir = param["folder"];
	namedir.append("round_");
	namedir.append("gamma_");
	namedir.append(douts(std::ceil(focals.gamma[0]*100)/100));
	namedir.append("regRew_");
	namedir.append(douts(std::ceil(focals.negReward[0] * 100) / 100)); 
	namedir.append("scaC_");
	namedir.append(douts(std::ceil(focals.scaleConst * 100) / 100));
	namedir.append("alphA_");
	namedir.append(douts(std::ceil(focals.alphaA * 100) / 100));
	namedir.append("alphC_");
	namedir.append(douts(std::ceil(focals.alphaC * 100) / 100));
	namedir.append("_seed");
	namedir.append(itos(param["seed"]));
	namedir.append(".txt");
	RoundOutput.open(namedir.c_str());
	if (param["data"] == "loc") {
		RoundOutput << "site_year	"
			<< "rel.abund.cleaners	"
			<< "rel.abund.visitors	"
			<< "rel.abund.residents	"
			<< "prob.Vis.Leav	"
			<< "market_binomial_data	"
			<< "market_binomial_pred	"
			<< "competence" << endl;
	}
	else{
		RoundOutput << "site_year	"
			<< "CleanerID	"
			<< "rel.abund.cleaners	"
			<< "rel.abund.visitors	"
			<< "rel.abund.residents	"
			<< "prob.Vis.Leav	"
			<< "visitorChoices	"
			<< "visitorChoices_pred	"
			<< "competence" << endl;
	}
	
}

void printRoundFile(ofstream& roundOut, const std::vector<locat_point>& emp_data) {
	for (int id_data_point = 0; id_data_point < emp_data.size(); ++id_data_point) {
		roundOut << emp_data[id_data_point].site_year << '\t' <<
			emp_data[id_data_point].rel_abund_clean << '\t' <<
			emp_data[id_data_point].rel_abund_visitors << '\t' <<
			emp_data[id_data_point].rel_abund_resid << '\t' <<
			emp_data[id_data_point].prob_Vis_Leav << '\t' <<
			emp_data[id_data_point].market_exp_success << '\t' <<
			emp_data[id_data_point].marketPred << endl;
	}
}

void printRoundFile(ofstream& roundOut, const std::vector<cleaner_point>& emp_data) {
	for (int id_data_point = 0; id_data_point < emp_data.size(); ++id_data_point) {
		roundOut << emp_data[id_data_point].site_year << '\t' <<
			emp_data[id_data_point].cleanerID << '\t' <<
			emp_data[id_data_point].rel_abund_clean << '\t' <<
			emp_data[id_data_point].rel_abund_visitors << '\t' <<
			emp_data[id_data_point].rel_abund_resid << '\t' <<
			emp_data[id_data_point].prob_Vis_Leav << '\t' <<
			emp_data[id_data_point].countVisitor << '\t' <<
			emp_data[id_data_point].marketPred << '\t' << 
			emp_data[id_data_point].group << endl;
	}
}

//double	getkforGamma(double& mode, double var) {
//	double k1 = (2 + (pow(mode, 2)) / var + 
//		mode * sqrt((4 + pow(mode, 2) / var) / var)) / 2;
//	double k2 = (2 + pow(mode, 2) / var - 
//		mode * sqrt((4 + pow(mode, 2) / var) / var)) / 2;
//	if (k1 > 0 || k2 < 0) return(k1);
//	else
//		if (k1 < 0 || k2 >0) return(k2);
//		else error("result out of range", CURRENT_FUNCTION);
//}


//double calcTransProbRatio(model_param focalParam, model_param newParam,
//	json &sim_param) {
//	double Pold2new = 1, Pnew2old=1, alphaBeta, 
//		betaBeta, kgamma,thetaGamma, ratio;
//	//gamma parameter
//	alphaBeta = getAlphaforBeta(newParam.gamma, sim_param["sdPert"][2]);
//	betaBeta = getbetaforBeta(newParam.gamma, sim_param["sdPert"][2]);
//	Pnew2old *= rnd::beta_pdf(focalParam.gamma, alphaBeta, betaBeta);
//	alphaBeta = getAlphaforBeta(focalParam.gamma, sim_param["sdPert"][2]);
//	betaBeta = getbetaforBeta(focalParam.gamma, sim_param["sdPert"][2]);
//	Pold2new *= rnd::beta_pdf(newParam.gamma, alphaBeta, betaBeta);
//	//NegReward parameter
//	kgamma = getkforGamma(newParam.negReward, sim_param["sdPert"][3]);
//	thetaGamma = getthetaforGamma(newParam.negReward, sim_param["sdPert"][3]);
//	Pnew2old *= rnd::pdf_gamma(focalParam.negReward, kgamma, thetaGamma);
//	kgamma = getkforGamma(focalParam.negReward, sim_param["sdPert"][3]);
//	thetaGamma = getthetaforGamma(focalParam.negReward, sim_param["sdPert"][3]);
//	Pold2new *= rnd::pdf_gamma(newParam.negReward, kgamma, thetaGamma);
//
//	ratio = Pnew2old / Pold2new;
//	return ratio;
//}


//double calculate_fit(const std::vector< data_point>& empData,
//	const std::vector< data_point>& simData) {
//	double sum_log_likelihood = 0.0;
//	for (int i = 0; i < empData.size(); ++i) {
//		sum_log_likelihood += log(simData[i].market_exp_success*empData[i].market_exp_success +
//			(1 - simData[i].market_exp_success)*(1 - empData[i].market_exp_success));
//	}
//	return(sum_log_likelihood);
//}

double calculate_fit(const std::vector<locat_point>& empData) {
	double sum_log_likelihood = 0.0;
	for (int i = 0; i < empData.size(); ++i) {
		sum_log_likelihood +=
			log(rnd::pdf_binomial(empData[i].countMarket, 
				empData[i].totalMarket, empData[i].marketPred));
	}
	return(sum_log_likelihood);
}

double calculate_fit(const std::vector<cleaner_point>& empData) {
	double sum_log_likelihood = 0.0;
	for (int i = 0; i < empData.size(); ++i) {
		sum_log_likelihood +=log(rnd::pdf_binomial(empData[i].countVisitor,
			empData[i].totalMarket, empData[i].marketPred))	;
	}
	return(sum_log_likelihood);
}

//model_param perturb_parameters(model_param focal_param, json &sim_param) {
//	model_param new_param;
//	// also, you can throw in your own random number generator that you want,
//	// I just add a number N(0, sd);
//	if (sim_param["pertScen"] == 0) {
//		new_param.alphaA = focal_param.alphaA + rnd::normal(0,
//			float(sim_param["sdPert"][0]));
//		new_param.alphaC = focal_param.alphaC + rnd::normal(0,
//			float(sim_param["sdPert"][1]));
//		new_param.negReward = focal_param.negReward;
//		new_param.gamma = focal_param.gamma;
//	}
//	if (sim_param["pertScen"] < 2) {
//		double alphaBeta = getAlphaforBeta(focal_param.gamma, sim_param["sdPert"][2]);
//		double betaBeta = getbetaforBeta(focal_param.gamma, sim_param["sdPert"][2]);
//		new_param.gamma = rnd::beta(alphaBeta, betaBeta);
//		double kgamma = getkforGamma(focal_param.negReward, sim_param["sdPert"][3]);
//		double thetaGamma = getthetaforGamma(focal_param.negReward,
//			sim_param["sdPert"][3]);
//		new_param.negReward = rnd::gamma(kgamma, thetaGamma);
//		new_param.alphaA = focal_param.alphaA;
//		new_param.alphaC = focal_param.alphaC;
//
//	}
//	else if (sim_param["pertScen"] == 2)
//	{
//		new_param.alphaA = focal_param.alphaA;
//		new_param.alphaC = focal_param.alphaC;
//		double alphaBeta = getAlphaforBeta(focal_param.gamma, sim_param["sdPert"][2]);
//		double betaBeta = getbetaforBeta(focal_param.gamma, sim_param["sdPert"][2]);
//		new_param.gamma = rnd::beta(alphaBeta, betaBeta);
//		new_param.negReward = focal_param.negReward;
//	}
//	else
//	{
//		new_param.alphaA = focal_param.alphaA;
//		new_param.alphaC = focal_param.alphaC;
//		double kgamma = getkforGamma(focal_param.negReward, sim_param["sdPert"][3]);
//		double thetaGamma = getthetaforGamma(focal_param.negReward,
//			sim_param["sdPert"][3]);
//		new_param.negReward = rnd::gamma(kgamma, thetaGamma);
//		new_param.gamma = focal_param.gamma;
//	}
//	return(new_param);
//}

double boundedParUnifPert(double & parVal, float pertRang ,double min, double max) {
	if (parVal - pertRang * 0.5 < min) {
		return (pertRang*rnd::uniform());
	}
	else if (parVal + pertRang * 0.5 > max){
		return (pertRang*rnd::uniform() + 1-pertRang);
	}
	else	{
		return (pertRang*(rnd::uniform() - 0.5) + parVal);
	}
}

model_param perturb_parameters_uniform(model_param focal_param, json &sim_param) {
	model_param new_param;
	// also, you can throw in your own random number generator that you want,
	// I just add a number N(0, sd);

	new_param.alphaA = bool(sim_param["pertScen"][0])*boundedParUnifPert(focal_param.alphaA,
		float(sim_param["sdPert"][0]), 0, INFINITY) +
		(!bool(sim_param["pertScen"][0]))*focal_param.alphaA;
	new_param.alphaC = bool(sim_param["pertScen"][0]) * boundedParUnifPert(focal_param.alphaC,
		float(sim_param["sdPert"][1]), 0, INFINITY) +
		(!bool(sim_param["pertScen"][1]))*focal_param.alphaC;
	new_param.gamma[0] = bool(sim_param["pertScen"][2]) *
		boundedParUnifPert(focal_param.gamma[0],
			sim_param["sdPert"][2], -1, 1) +
			(!bool(sim_param["pertScen"][2]))*focal_param.gamma[0];
	new_param.negReward[0] = bool(sim_param["pertScen"][3])*
		boundedParUnifPert(focal_param.negReward[0],
			sim_param["sdPert"][3], -INFINITY, INFINITY) +
			(!bool(sim_param["pertScen"][3]))*focal_param.negReward[0];
	if (sim_param["Group"]) {
		new_param.gamma[1] = bool(sim_param["pertScen"][2]) *
			boundedParUnifPert(focal_param.gamma[1],
				sim_param["sdPert"][2], -1, 1) +
				(!bool(sim_param["pertScen"][2]))*focal_param.gamma[1];
		new_param.negReward[1] = bool(sim_param["pertScen"][3])*
			boundedParUnifPert(focal_param.negReward[1],
				sim_param["sdPert"][3], -INFINITY, INFINITY) +
				(!bool(sim_param["pertScen"][3]))*focal_param.negReward[1];
	}
	new_param.scaleConst = bool(sim_param["pertScen"][4])* boundedParUnifPert(focal_param.scaleConst,
		sim_param["sdPert"][4], 0, INFINITY) +
		(!bool(sim_param["pertScen"][4]))*focal_param.scaleConst;
		/*bool(sim_param["pertScen"][4])* rnd::normal(0, sim_param["sdPert"][4])+
		focal_param.scaleConst;*/
	return(new_param);
}

void do_simulation(//del focal_model,
	std::vector<locat_point> &emp_data, model_param focal_comb,
	json sim_param) {
	client *clientSet;
	clientSet = new client[int(sim_param["totRounds"]) * 2];
	int idClientSet;
	FAATyp1 Cleaner (focal_comb.alphaC, focal_comb.gamma[0],
		focal_comb.negReward[0],focal_comb.alphaA);
	double VisPref, init;
	int countRVopt;
	// Loop through the data points
	for (int id_data_point = 0; id_data_point < emp_data.size(); ++id_data_point) {
		init = focal_comb.gamma[0]*
			(1 - pow(1 -
				emp_data[id_data_point].rel_abund_resid -
				emp_data[id_data_point].rel_abund_visitors, 2)) / (1 - focal_comb.gamma[0]);
		Cleaner.rebirth(init);
		draw(clientSet, sim_param["totRounds"],
			emp_data[id_data_point].rel_abund_resid,
			emp_data[id_data_point].rel_abund_visitors);
		idClientSet = 0;
		VisPref = 0, countRVopt = 0;
		// Loop through the learning rounds
		for (int trial = 0; trial < sim_param["totRounds"]; ++trial) {
			Cleaner.act(clientSet, idClientSet, emp_data[id_data_point].prob_Vis_Leav,
				sim_param["ResProbLeav"], sim_param["VisReward"],
				sim_param["ResReward"], sim_param["inbr"], sim_param["outbr"],
				learnScenario(sim_param["scenario"]));
			Cleaner.update();
			if (trial > int(sim_param["totRounds"]) * float(sim_param["propfullPrint"])) {
				if (Cleaner.getstate(0) == 0) {
					++countRVopt;
					if (Cleaner.cleanOptionsT[Cleaner.getChoice(0)] == visitor) ++VisPref;
				}
			}
		}
		if (countRVopt == 0) emp_data[id_data_point].marketPred = 0.5;
		else emp_data[id_data_point].marketPred = VisPref / countRVopt;
		Cleaner.rebirth();
		}
	delete[] clientSet;
	return;
}


void do_simulation(//del focal_model,
	std::vector<cleaner_point> &emp_data, model_param focal_comb,
	json sim_param) {
	client *clientSet;
	clientSet = new client[int(sim_param["totRounds"]) * 2];
	int idClientSet;
	FAATyp1 *cleaners[2];
	cleaners[0] = new FAATyp1(focal_comb.alphaC, focal_comb.gamma[0],
		focal_comb.negReward[0], focal_comb.alphaA);
	if(sim_param["Group"]) cleaners[1] = new FAATyp1(focal_comb.alphaC, focal_comb.gamma[1],
		focal_comb.negReward[1], focal_comb.alphaA);
	double VisPref, init;
	int countRVopt;
	abs2rel_abund(emp_data, focal_comb);
	// Loop through the data points
	for (int id_data_point = 0; id_data_point < emp_data.size(); ++id_data_point) {
		if (id_data_point > 0 &&
			(emp_data[id_data_point].site_year == emp_data[id_data_point - 1].site_year &&
			emp_data[id_data_point].group == emp_data[id_data_point - 1].group)) {
			emp_data[id_data_point].marketPred = emp_data[id_data_point - 1].marketPred;
		}
		else
		{
			init = focal_comb.gamma[emp_data[id_data_point].group]*
				(1 - pow(1 -
					emp_data[id_data_point].rel_abund_resid -
					emp_data[id_data_point].rel_abund_visitors, 2)) / 
					(1 - focal_comb.gamma[emp_data[id_data_point].group]);
			cleaners[emp_data[id_data_point].group]->rebirth(init);
			draw(clientSet, sim_param["totRounds"],
				emp_data[id_data_point].rel_abund_resid,
				emp_data[id_data_point].rel_abund_visitors);
			idClientSet = 0;
			VisPref = 0, countRVopt = 0;
			// Loop through the learning rounds
			for (int trial = 0; trial < sim_param["totRounds"]; ++trial) {
				cleaners[emp_data[id_data_point].group]->
					act(clientSet, idClientSet, emp_data[id_data_point].prob_Vis_Leav,
					sim_param["ResProbLeav"], sim_param["VisReward"],
					sim_param["ResReward"], sim_param["inbr"], sim_param["outbr"],
					learnScenario(sim_param["scenario"]));
			cleaners[emp_data[id_data_point].group]->update();
				if (trial > int(sim_param["totRounds"]) * float(sim_param["propfullPrint"])) {
					if (cleaners[emp_data[id_data_point].group]->getstate(0) == 0) {
						++countRVopt;
						if (cleaners[emp_data[id_data_point].group]->
							cleanOptionsT[cleaners[emp_data[id_data_point].group]->getChoice(0)] == visitor) 
							++VisPref;
					}
				}
			}
			if (countRVopt == 0) emp_data[id_data_point].marketPred = 0.5;
			else emp_data[id_data_point].marketPred = VisPref / countRVopt;
			cleaners[emp_data[id_data_point].group]->rebirth();
		}
	}
	delete[] clientSet;
	return;
}

void checkGroups(vector<cleaner_point> &data,bool groups) {
	if (!groups) {
		vector<cleaner_point>::iterator dataIt;
		for (dataIt = data.begin(); dataIt < data.end(); ++dataIt)
			dataIt->group = 0;
	}
	return;
}

int main(int argc, char* argv[]){

	mark_time(1);

	// Hardwire parameter values:
	// Only for debugging 
	// input parameters provided by a JSON file with the following
	// structure:
	//json sim_param;
	//sim_param["totRounds"]    = 5000;
	//sim_param["ResReward"]    = 1;
	//sim_param["VisReward"]    = 1;
	//sim_param["ResProbLeav"]  = 0;
	//sim_param["scenario"]  = 0;
	//sim_param["inbr"]         = 0;
	//sim_param["outbr"]        = 0;
	//sim_param["seed"]         = 3;
	//sim_param["forRat"]       = 0.0;
	//sim_param["propfullPrint"]       = 0.7;
	//sim_param["sdPert"]       = {0.05, 0.05 ,0.15 ,0.1, 10}; 
	//// alphaA, alphaC, Gamma, NegRew,scaleConst
	//sim_param["chain_length"] = 100;
	//sim_param["init"]       = {0.05, 0.05 , 0.93,0.02, 58};
	//sim_param["init2"] =	{ 0.05, 0.05 , 0.93,0.02, 58 };
	// //alphaA, alphaC, gamma, NegRew, scaleConst
	//sim_param["pertScen"] = {false,false,true,true,true};
	////enum perturnScen {all,  bothFut, justGam, justNegRew};
	//sim_param["MCMC"] = 0;
	//sim_param["data"] = "clean"; // "loc", "clean"
	//sim_param["nRep"] = 30 ;
	//sim_param["folder"] = "M:/Projects/LearnDataModel/Simulations/test_/";
	//sim_param["dataFile"] = "M:/Projects/LearnDataModel/Data/data_cleaner_abs_threa1.5.txt";
	//sim_param["Group"] = true;

	////ifstream marketData ("E:/Projects/Clean.ActCrit/Data/data_ABC.txt");
	

	// reading of parameters: 
	ifstream parameters(argv[1]);
	if (parameters.fail()) { cout << "JSON file failed" << endl; }
	json sim_param = nlohmann::json::parse(parameters);
	
	// Set random seed
	rnd::set_seed(sim_param["seed"]);
	
	//vector< data_point > emp_data = read_Data(marketData); 
	ifstream marketData(sim_param["dataFile"].get<std::string>());//("I:/Projects/Clean.ActCrit/Data/data_ABC_site.txt");
	//marketData.open(, ifstream::in);
	vector <locat_point> emp_data_loc;
	//ifstream marketData_cleaner; // ("I:/Projects/Clean.ActCrit/Data/data_ABC_cleaner.txt");
	vector <cleaner_point> emp_data_clean; 
	
	if (sim_param["data"] == "loc") emp_data_loc = read_locData(marketData);
	else emp_data_clean = read_CleanData(marketData);

	checkGroups(emp_data_clean, sim_param["Group"]);

	// read the data structured by locations
	model_param init_parameters; 
	init_parameters.alphaA = sim_param["init"][0];
	init_parameters.alphaC = sim_param["init"][1];
	init_parameters.gamma[0] = sim_param["init"][2];
	init_parameters.negReward[0] = sim_param["init"][3];
	init_parameters.gamma[1] = sim_param["init"][2];
	init_parameters.negReward[1] = sim_param["init"][3];
	init_parameters.scaleConst = sim_param["init"][4];

	model_param focal_param = init_parameters;
	
	if (sim_param["MCMC"] == 1) {

		double curr_loglike;

		// we calculate the fit of the starting point, first we simulate data using 
		// the initial parameters
		// we also pass on the empirical data, to use the x and y coordinates.
		if (sim_param["data"] == "loc") {
			do_simulation(//focal_model
				emp_data_loc, init_parameters, sim_param);
			// function that calculates fit
			curr_loglike = calculate_fit(emp_data_loc);
			while (isinf(-curr_loglike)) {
				focal_param = perturb_parameters_uniform(focal_param, sim_param);
				curr_loglike  = calculate_fit(emp_data_loc);
			}
		}
		else {
			do_simulation(//focal_model
			emp_data_clean, init_parameters, sim_param);
			curr_loglike = calculate_fit(emp_data_clean);
			while (isinf(-curr_loglike)) {
				focal_param = perturb_parameters_uniform(focal_param, sim_param);
				curr_loglike  = calculate_fit(emp_data_clean);
			}
		}


		ofstream outfile;
		initializeChainFile(outfile, sim_param);
		double new_loglike, ratio;
		for (int r = 0; r < sim_param["chain_length"]; ++r) {  // 
			//cout << "Iteration	" << r << endl;
			model_param new_param = perturb_parameters_uniform(focal_param, sim_param);
			//if ((new_param.gamma < 1 && new_param.gamma > -1) &&
			//	((new_param.negReward <= INFINITY && new_param.negReward >= -INFINITY) &&
			//	(new_param.scaleConst <= INFINITY && new_param.scaleConst > 0))) {
			if (sim_param["data"] == "loc") {
				do_simulation(//focal_model, 
					emp_data_loc, focal_param, sim_param);
				curr_loglike = calculate_fit(emp_data_loc);
				do_simulation(//focal_model, 
					emp_data_loc, new_param, sim_param);
				new_loglike = calculate_fit(emp_data_loc);
			}
			else {
				do_simulation(//focal_model, 
					emp_data_clean, focal_param, sim_param);
				curr_loglike = calculate_fit(emp_data_clean);
				do_simulation(//focal_model, 
					emp_data_clean, new_param, sim_param);
				new_loglike = calculate_fit(emp_data_clean);
			}
				//ratio = 1;// calcTransProbRatio(focal_param, new_param, sim_param);
			/*}
			else
			{
				new_loglike = -INFINITY, ratio = 1;
			}*/ 
			ratio = exp(new_loglike - curr_loglike);
			if(ratio==INFINITY) 
				wait_for_return();
			// better fit is larger, so ratio is > 1, so accept all.
			if (rnd::uniform() < ratio) {
				focal_param = new_param;
				curr_loglike = new_loglike;
			}
			outfile << r << "\t";
			outfile << focal_param.alphaA << "\t"
				<< focal_param.alphaC << "\t"
				<< focal_param.gamma[0] << "\t"
				<< focal_param.negReward[0] << "\t";
			if(sim_param["Group"])
				outfile << focal_param.gamma[1] << "\t"
				<< focal_param.negReward[1] << "\t";
			outfile	<< focal_param.scaleConst << "\t"
				<< curr_loglike << "\t";
			outfile << ratio << endl;
		}
		outfile.close();
		// done!
		//wait_for_return();
	}
	else {
		if (sim_param["Group"]) {
			init_parameters.gamma[1] = sim_param["init2"][2];
			init_parameters.negReward[1] = sim_param["init2"][3];
		}
		if (sim_param["data"] == "loc") {
			ofstream roundOut;
			initializeRoundFile(roundOut, sim_param, init_parameters);
			double average_prediction[12];
			for (int j = 0; j < 12; ++j) average_prediction[j] = 0;
			for (int rep = 0; rep < int(sim_param["nRep"]); ++rep) {
				do_simulation(emp_data_loc, init_parameters, sim_param);
				for (int foc_datapoint = 0; foc_datapoint < emp_data_loc.size();
					++foc_datapoint) {
					average_prediction[foc_datapoint] +=
						emp_data_loc[foc_datapoint].marketPred;
				}
			}
			for (int j = 0; j < emp_data_loc.size(); ++j)
				emp_data_loc[j].marketPred = 
				  average_prediction[j] / double(sim_param["nRep"]);
			printRoundFile(roundOut, emp_data_loc);
			roundOut.close();
		}
		else {
			ofstream roundOut;
			initializeRoundFile(roundOut, sim_param, init_parameters);
			double average_prediction[120];
			for (int j = 0; j < 120; ++j) average_prediction[j] = 0;
			for (int rep = 0; rep < int(sim_param["nRep"]); ++rep) {
				do_simulation(emp_data_clean, init_parameters, sim_param);
				for (int foc_datapoint = 0; foc_datapoint < emp_data_clean.size();
					++foc_datapoint) {
					average_prediction[foc_datapoint] +=
						emp_data_clean[foc_datapoint].marketPred;
				}
			}
			for (int j = 0; j < emp_data_clean.size(); ++j) 
				emp_data_clean[j].marketPred = 
				  average_prediction[j]/double(sim_param["nRep"]);
			printRoundFile(roundOut, emp_data_clean);
			roundOut.close();
		}
		
	}
	mark_time(0);
	//wait_for_return();
	return 0;
}