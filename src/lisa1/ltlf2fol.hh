/*-------------------------------------------------------------------*/
// Yong Li (liyong@ios.ac.cn)
/*-------------------------------------------------------------------*/

#pragma once

#include <sstream>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "spotutil.hh"



/*-------------------------------------------------------------------*/
//Translate an LTLf formula to a First-Order Logic (FOL) formula, as described in 
// the paper https://www.cs.rice.edu/~vardi/papers/ijcai13.pdf
/*-------------------------------------------------------------------*/
void
ltlf_to_fol(ostream& os, formula& f);
// helper for translation from LTLf formulas to FOL formulas
string
translate2fol(formula& f, int t, int& c);

/*-------------------------------------------------------------------*/
//Translate an LTLf formula to a past First-Order Logic (FOL) formula, as described in 
// the paper https://arxiv.org/abs/1901.06108 by Shufang Zhu et al.
/*-------------------------------------------------------------------*/
void
ltlf_to_pfol(ostream& os, formula& f);
// helper for translation from LTLf formulas to FOL formulas
string
translate2pfol(formula& f, string t, int& c);


// TODO: no idea about what this function will do 
void
trans_prefixltlf2fol(ostream &os, formula &f);

string
get_prefix2fol(formula& f, int t, int& c);



/*-------------------------------------------------------------------*/
// Normal form
/*-------------------------------------------------------------------*/

// obtain a formula in negation normal form (NNF)
formula
get_nnf(formula& f);

// helper for obataining NNF for negation
formula
push_not_in(formula& f);

// obtain a formula in Boolean normal form (BNF)
// BNF is choosed by default
formula
get_bnf(formula& f);



