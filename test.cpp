#include <Rcpp.h>
#include "../Cpp/json.hpp"  
using json = nlohmann::json;
using namespace Rcpp;


// [[Rcpp::depends(RcppJson)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

std::string hello() {
  return "hello";
}
void bla() {
  Rprintf("hello\\n");
}
void bla2( int x, double y) {
  Rprintf("hello (x = %d, y = %5.2f)\\n", x, y);
}
class World {
public:
  World() : msg("hello") {}
  void set(std::string msg) { this->msg = msg; }
  std::string greet() { return msg; }
private:
  std::string msg;
};

RCPP_MODULE(yada){
  using namespace Rcpp;
  function("hello" , &hello);
  function("bla" , &bla);
  function("bla2" , &bla2);
  class_<World>("World")
    .constructor()
    .method("greet", &World::greet)
    .method("set", &World::set)
  ;
}



struct model_param {
  model_param()=default;
  double logist() {
    return (1 / (1 + exp(-(alphaA-alphaC))));
  }
  void set(double alphaC_,double alphaA_){
    this->alphaA = alphaA_;
    this->alphaC = alphaC_;
  }
  void copy(Rcpp::XPtr<model_param> mp_ptr){
    this->alphaA = mp_ptr->alphaA;
    this->alphaC = mp_ptr->alphaC;
  }
  Rcpp::XPtr<model_param>get_ptr(){
    Rcpp::XPtr<model_param> p(this,true);
    return p;
  }
  double alphaC, alphaA;
  json myJson = json::parse(R"({
    "pi": 3.141,
    "happy": true
  }
  )");
};



// [[Rcpp::export]]
double showPi(json myJson){
  return myJson["pi"];
}


RCPP_MODULE(cleaner){
  using namespace Rcpp;
  function("timesTwo", &timesTwo);

  class_<model_param>("model_param")
    .constructor()
    .method("set_gamma", &model_param::set)
    .method("logist", &model_param::logist)
    .method("copy", &model_param::copy)
    .method("get_ptr",&model_param::get_ptr)
  ;
}



// struct model_param {
//   //model_param(model_param const &obj);
//   model_param()=default;
//   double alphaC, alphaA, scaleConst;
//   double gamma[2], negReward[2];
//   double probFAA[2];
//   double interpReg, slopRegRelAC, slopRegPVL;
//   void copy (Rcpp::XPtr<model_param> mp_ptr) {
//     alphaA = mp_ptr->alphaA;
//     alphaC = mp_ptr->alphaC;
//     scaleConst = mp_ptr->scaleConst;
//     gamma[0] = mp_ptr->gamma[0];
//     negReward[0] = mp_ptr->negReward[0];
//     probFAA[0] = mp_ptr->probFAA[0];
//     gamma[1] = mp_ptr->gamma[1];
//     negReward[1] = mp_ptr->negReward[1];
//     probFAA[1] = mp_ptr->probFAA[1];
//   }
//   Rcpp::XPtr<model_param>get_ptr(){
//     Rcpp::XPtr<model_param> p(this,true);
//     return p;
//   }
// };
// 
// RCPP_MODULE(cleaner){
//   using namespace Rcpp;
// 
//   class_<model_param>("model_param")
//     .constructor()
//     .method("copy", &model_param::copy)
//     .method("get_ptr",&model_param::get_ptr)
//   ;
// }

class Uniform {
public:
  Uniform(double min_, double max_) :
  min(min_), max(max_) {}
  NumericVector draw(int n) {
    RNGScope scope;
    return runif(n, min, max);
  }
private:
  double min, max;
};

RCPP_MODULE(UNIF){
  using namespace Rcpp;
  class_<Uniform>("Uniform")
    .constructor<double,double>()
    .method("draw", &Uniform::draw)
  ;
  
}