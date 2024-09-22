template<class Type>
Type erf(Type x){
  Type y = sqrt(Type(2.0)) * x;
  Type z = Type(2.0) * pnorm(y) - Type(1.0);
  return z;
}
