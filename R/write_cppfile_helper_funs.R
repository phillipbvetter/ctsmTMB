########## GENERAL HELPER FUNCTIONS ###############
############################################################
write_helper_cppfunctions = function(){
  
  txt = c()
  
  ##################################################
  # Find function finds indices of NAs in a vector
  ##################################################
  newtxt = "\n//////////// FIND NA INDICES IN VECTOR ///////////
  template <class Type>
  vector<Type> is_not_na(vector<Type> x){
    vector<Type> y(x.size());
    y.fill(Type(1.0));
      for(int i=0; i<x.size(); i++){
        if( R_IsNA(asDouble(x(i))) ){
          y(i) = Type(0.0);
        }
      }
    return y;
  }"
  txt = c(txt,newtxt)
  
  ##################################################
  # This function removes NAs from a vector
  ##################################################
  newtxt = "\n//////////// REMOVE NA'S FROM VECTOR ///////////
  template<class Type>
  vector<Type> remove_nas__(vector<Type> obsVec, int number_of_nonNA_observations, vector<Type> is_not_na_vector){
    int ii = 0;
    vector<Type> y_reduced(number_of_nonNA_observations);
      for(int i=0; i < obsVec.size(); i++){
        if(is_not_na_vector(i) == Type(1.0)){
          y_reduced(ii) = obsVec(i);
          ii++;
        }
      }
    return y_reduced;
  }"
  txt = c(txt,newtxt)
  
  ##################################################
  # Construct Permutation Matrix E for Kalman Filter NA-removal
  ##################################################
  newtxt = "\n//////////// helper fun: construct permutation matrix ///////////
  template <class Type>
  matrix<Type> construct_permutation_matrix(int number_of_nonNA_observations, int number_of_obs_eqs, vector<Type> is_not_na_vector){
	  matrix<Type> E(number_of_nonNA_observations, number_of_obs_eqs);
	  E.setZero();
	  int j=0;
	  for(int i=0; i < number_of_obs_eqs; i++){
      /*if p(i) is 1 then include by setting 1 in diagonal of matrix*/
		  if(is_not_na_vector(i) == Type(1.0)){
        E(j,i) = Type(1.0);
			  j += 1;
      }
	  }
	  return E;
  }"
  txt = c(txt,newtxt)
  
  ##################################################
  # Implements Tukey and Huber loss functions
  ##################################################
  newtxt = "\n//////////// loss function ///////////
  template<class Type>
  Type lossfunction__(Type x, vector<Type> tukeypars, Type huber_c, int lossFunc){
    Type loss;
    if(lossFunc==1){
      Type a = tukeypars(0);
      Type b = tukeypars(1);
      Type c = tukeypars(2);
      Type d = tukeypars(3);
      loss = d * ( (Type(1.0)/(Type(1.0)+exp(-a*(x-b)))) + c );
    } else if (lossFunc==2){
      Type c_squared = pow(huber_c,2);
      loss = c_squared * (sqrt(1 + (x / c_squared)) - 1);
    } else {
      loss = x;
    }
    return loss;
  }"
  txt = c(txt,newtxt)
  
  ##################################################
  # MAP estimation
  ##################################################
  newtxt = "\n//////////// MAP estimation helper ///////////
  template<class Type>
  vector<Type> get_free_pars__(vector<int> mapints, int sum_mapints, vector<Type> parvec) {
	  vector<Type> ans(sum_mapints);
	  int j=0;
	  for(int i=0;i<mapints.size();i++){
		  if(mapints(i)==1){
			  ans(j) = parvec(i);
			  j += 1;
		  }
	  }
	  return ans;
  }"
  txt = c(txt, newtxt)
  
  return(txt)
}

########## CUSTOM FUNCTIONS ###############
############################################################
write_custom_functions = function(){
  
  txt = c()
  
  newtxt = "\n//////////// Error Function ///////////
  template<class Type>
  Type erf(Type x){
    Type y = sqrt(Type(2.0)) * x;
    Type z = Type(2.0) * pnorm(y) - Type(1.0);
    return z;
    }"
  txt = c(txt, newtxt)
  
  return(txt)
}
