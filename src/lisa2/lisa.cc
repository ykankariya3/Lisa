#include <algorithm>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <unordered_map>
#include <string>
#include <vector>
#include <queue>
#include <list>
#include <set>
#include <chrono> 

// spot 
#include <spot/twaalgos/hoa.hh>
#include <spot/twaalgos/translate.hh>
#include <spot/twaalgos/isdet.hh>
#include <spot/twaalgos/product.hh>
#include <spot/twaalgos/complement.hh>
#include <spot/twaalgos/sum.hh>
#include <spot/twaalgos/word.hh>
#include <spot/twaalgos/degen.hh>
#include <spot/twaalgos/remprop.hh>
#include <spot/twaalgos/minimize.hh>

#include <spot/twa/twagraph.hh>

#include <spot/tl/formula.hh>
#include <spot/tl/print.hh>
#include <spot/tl/parse.hh>
#include <spot/tl/relabel.hh>
#include <spot/tl/ltlf.hh>
#include <spot/tl/simplify.hh>

#include <spot/misc/optionmap.hh>
#include <spot/misc/timer.hh>

#include "mona.hh"
#include "spotutil.hh"
#include "dfwavar.hh"
#include "dfwa.hh"
#include "debug.hh"
#include "synt.hh"
#include "dfwamin.hh"
#include "dfwamin2.hh"
//#include "dfwamin3.hh"

using namespace spot;
using namespace std;

#define info(o) cout << "[INFO] " << ( o ) << endl
#define erro(o) cerr << "[ERRO] " << ( o ) << endl

// options

static struct opt_t
{
	const char* _ltlfile_name = nullptr;
	const char* _parfile_name = nullptr;

	bool _symbolic = true;
	bool _minimization = false;

	unsigned _num_ap_for_mona = 7;
	unsigned _num_product = 6;
	unsigned _num_st_for_single = 800;
	unsigned _num_st_for_product = 2500;
	int _num_last_automata = -1;

	bool _synthesis = false;
	bool _out_start = false;
	bool _env_first = false;

	uint8_t _bdd = 0;

}* opt;


void 
get_formulas(vector<formula>& lst, formula f)
{
    cout << "Breaking formula into small pieces..." << endl;
    if(f.kind() == op::And)
    {
        // needs to limit the number of conjunctions if no minimization is used before producting two FAs
        for(formula child: f)
          {
              lst.push_back(child);
              //cout << "subformula: " << child << endl;
          }
    }else
    {
        lst.push_back(f);
    }
    //cout << "split formula " << endl;
}

class dfwa_pair
{
public:
	unsigned _num_states;
	bool _is_explicit;
	twa_graph_ptr _twa;
	dfwa* _dfa = nullptr;
	unsigned _num_propduct = 0;

	formula _formula;

	dfwa_pair(twa_graph_ptr aut, unsigned num_states, bool is_explicit, formula& f)
	: _num_states(num_states), _is_explicit(is_explicit), _formula(f)
	{
		_twa = aut;
	}
	dfwa_pair(dfwa* aut, unsigned num_states, bool is_explicit, formula& f)
		: _num_states(num_states), _is_explicit(is_explicit), _formula(f)
	{
		_dfa = aut;
		//_twa = new shared_ptr<twa_graph>(new twa_graph);
	}
};

struct GreaterThanByDfwaSize
{
  bool operator()(dfwa_pair& p1, dfwa_pair& p2) const
  {
    if(p1._num_states < p2._num_states)
    {
        return false;
    }else
    if(p1._num_states == p2._num_states)
    {
    	return !p1._is_explicit;
    }
    return p1._num_states >= p2._num_states;
  }
};

bool 
compare_aut_size(twa_graph_ptr p1, twa_graph_ptr p2)
{
    if(p1->num_states() == p2->num_states())
    {
        return false;
    }
    return p1->num_states() > p2->num_states();
}

twa_graph_ptr minimize_explicit(twa_graph_ptr A)
{
	twa_graph_ptr C = spot::minimize_wdba(A);
	//A = spot::minimize_obligation(A);
	// check equivalence of two automata
#ifdef DEBUG
	string word = is_twa_equivalent(A, C);
	if(word.size() == 0)
	{
		cout << "A: equivalent two automata" << endl;
	}
#endif
	return C;
}

tuple<dfwa*, unsigned, bool>
minimize_symbolic(dfwa_ptr P)
{
	unsigned num_states;
	tuple<dfwa*, unsigned, bool> result;
	if(opt->_bdd == 1)
	{
		cudd manager;
		cudd_ptr mgr = &manager;
		// call dfa minimization
		dfwa_min min(mgr, P);
		min.minimize();
		dfwa_ptr res = min.move_dfwa();
		// check wehter the emptiness is changed
		num_states = min.get_num_min_states();
		result = make_tuple<>(&res, num_states, true);
	}else
	if(opt->_bdd == 2)
	{
		dfwa_min_bdd min(P);
		min.minimize();
		dfwa_ptr res = min.move_dfwa();
		// check wehter the emptiness is changed
		num_states = min.get_num_min_states();
		result = make_tuple<>(&res, num_states, true);
	}else
	{
		/*
		dfwa_min_sylvan min(P);
		min.minimize();
		dfwa_ptr res = min.move_dfwa();
		// check wehter the emptiness is changed
		num_states = min.get_num_min_states();
		result = make_tuple<>(&res, num_states, true);
	
		*/
	}
	return result;
}

tuple<dfwa*, unsigned, bool>
make_product(bdd_dict_ptr dict, dfwa* A, dfwa* B, unsigned num_prod)
{
	unsigned num_states;

	dfwa_ptr P = product_dfwa_and(*A, *B);
	//cout << "labels in product and: " << P._label_cube << endl;
	//cout << "state_0 in product and: " << P._curr_cube << endl;
	//cout << "state_1 in product and: " << P._next_cube << endl;
	//cout << "product: " << endl;
	//P.output(cout);

	//unsigned var_num = 1;
	//cout << "condition: " << (P._state_vars.get_var_num(0) > var_num) << endl;

	if(num_prod > (opt->_num_product))
	{
		tuple<dfwa*, unsigned, bool> result = minimize_symbolic(P);
		delete &P;
		cout << "return from minimal product..." << endl;
		return result;
	}else
	{
		num_states = bdd_nodecount(P._trans);
		return make_tuple<>(&P, num_states, false);
	}

}



void print_usage()
{
	cout << "Usage: lisa [OPTION...] [FILENAME[/COL]...]" << endl;
	cout << "Read a formula file and output the number of states of the constructed DFA" << endl << endl;
	cout << " Input options:" << endl;
	cout << " -h  " << "                  show this help page" << endl;
	cout << " -exp" << "                  use only explicit method (default false)" << endl;
	cout << " -min" << "                  minimize the last symbolic DFA (default false)" << endl;
	cout << " -syn" << "                  synthesize after DFA construction (default false)" << endl;
	cout << " -bdd" << "                  use buddy for DFA minimization" << endl;
	//cout << " -syl" << "                  use sylvan for DFA minimization (default)" << endl;
	cout << " -cdd" << "                  use cudd for DFA minimization" << endl;
	cout << " -nap" << "  <int>           number of atomic propositions for calling mona (default 7)" << endl;
	cout << " -npr" << "  <int>           number of products for calling minimization (default 6)" << endl;
	cout << " -nia" << "  <int>           number of states of individual DFA for calling symbolic approach (default 800)" << endl;
	cout << " -npa" << "  <int>           number of states of product DFA for calling symbolic approach (default 2500)" << endl;
	cout << " -lst" << "  <int>           number of last automata for calling symbolic approach (default -1)" << endl;
	cout << " -out" << "                  print out the wining strategy if realizable" << endl;
	cout << " -part" << " <file>          the file specifying the input and output propositions" << endl;
	cout << " -ltlf" << " <file>          the file specifying the input LTLf formula" << endl;
	cout << " -env" << "                  environment plays first" << endl;
}

void parse_opt(int argc, char** argv)
{
	// first one is lisa, separated by space
	if(argc == 1)
	{
		print_usage();
	}
	for(int i = 1; i < argc; i ++)
	{
		string s(argv[i]);
		//cout << argv[i] << endl;
		if(s.size() == 0)
		{
			continue;
		}
		if(s == "-exp")
		{
			opt->_symbolic = false;
			//cout << "hello" << s << endl;
			continue;
		}
		if(s == "-min")
		{
			opt->_minimization = true;
			continue;
		}
		if(s == "-syn")
		{
			opt->_synthesis = true;
			continue;
		}
		if(s == "-out")
		{
			opt->_out_start = true;
			continue;
		}
		if(s == "-env")
		{
			opt->_env_first = true;
			continue;
		}
		if(s == "-nap" && i + 1 < argc)
		{
			opt->_num_ap_for_mona = stoi(argv[i + 1]);
			i ++;
			continue;
		}
		if(s == "-npr" && i + 1 < argc)
		{
			opt->_num_product = stoi(argv[i + 1]);
			i ++;
			continue;
		}
		if(s == "-nia" && i + 1 < argc)
		{
			opt->_num_st_for_single = stoi(argv[i + 1]);
			i ++;
			continue;
		}
		if(s == "-lst" && i + 1 < argc)
		{
			opt->_num_last_automata = stoi(argv[i + 1]);
			i ++;
			continue;
		}
		if(s == "-npa" && i + 1 < argc)
		{
			opt->_num_st_for_product = stoi(argv[i + 1]);
			i ++;
			continue;
		}
		if(s == "-ltlf" && i + 1 < argc)
		{
			opt->_ltlfile_name = argv[i + 1];
			//cout << "hello" << argv[i+1] << endl;
			i ++;
			continue;
		}
		if(s == "-part" && i + 1 < argc)
		{
			opt->_parfile_name = argv[i + 1];
			i ++;
			continue;
		}
		if(s == "-cdd")
		{
			opt->_bdd = 1;
			continue;
		}
		if(s == "-bdd")
		{
			opt->_bdd = 2;
			continue;
		}
		/*
		if(s == "-syl")
		{
			opt->_bdd = 0;
			continue;
		}
		*/
		if(s == "-h")
		{
			print_usage();
			exit(0);
		}else
		{
			erro("wrong input options: " + s);
			print_usage();
			exit(-1);
		}
	}
	// validity checking
	if(opt->_ltlfile_name == nullptr )
	{
		erro( "missing LTLf file name");
		exit(-1);
	}
	if(opt->_synthesis && ( opt->_parfile_name == nullptr))
	{
		erro("missing proposition partition file name");
		exit(-1);
	}

}

dfwa*
symbolize_twa(bdd_dict_ptr dict, twa_graph_ptr aut)
{
	// there is alive states
	bdd label_cube = bddtrue;
	for (auto f : aut->ap())
	{
		bdd f_var = bdd_ithvar(aut->register_ap(f));
		label_cube = label_cube & f_var;
		//cout << "formula : " << f << " index: " << aut->register_ap(f) << endl;
	}

	set<unsigned> finals_aut;
	compute_final_states(aut, finals_aut);
	/*
	cout << "final states: " << endl;
	for(unsigned k : finals_aut)
	{
		cout << "final: " << k << endl;
	}*/
	/*
	if(aut->num_states() < 20)
	{
		ofstream outfile("output" + to_string(number) + ".hoa");
		print_hoa(outfile, aut);
		for(unsigned k : finals_aut)
		{
			cout << "output: " << k << endl;
		}
		number ++;
	}*/
	// now compute dfwa
	dfwa* A = new dfwa(aut, label_cube, finals_aut);
	return A;
}
void
read_from_part_file(const char *file_name, vector<string>& input, vector<string>& output)
{
    // const char * file_name
    ifstream part_file(file_name);
    if (part_file.is_open())
    {
        bool flag = false;
        string line;
        while(getline(part_file, line))
        {
            if(str_contain(line, "inputs"))
            {
                string delimiter = ":";
                line = line.substr(line.find(delimiter) + 1);
                str_split(line, input, ' ');
            }else
            if(str_contain(line, "outputs"))
            {
                string delimiter = ":";
                line = line.substr(line.find(delimiter) + 1);
                str_split(line, output, ' ');
            }else
            {
                cout << "read partfile error!" <<endl;
                cout << file_name <<endl;
                cout << line <<endl;
                exit(-1);
            }
        }
    }
}

std::string opToString(spot::op op)
{
    switch (op)
    {
        case spot::op::Not: return "!";
        case spot::op::And: return "&&";
        case spot::op::Or: return "||";
        default: return "";
    }
}

spot::bdd_dict_ptr dict = spot::make_bdd_dict();
int sizes = 0;

struct TreeNode {
    spot::op value;
    //spot::formula leafValue;
    std::vector<TreeNode*> children;
	dfwa_pair* dfa;

    TreeNode(const spot::op& val) : value(val) {}
};

TreeNode* createTree(spot::formula& f)
{
    TreeNode* node = new TreeNode(f.kind());
    
    if (node->value == op::And || node->value == op::Or)
    {
        for (formula children: f)
        {
            TreeNode* child = createTree(children);
            node->children.push_back(child);
        }
    } 
	// else if (node-> value == op::Not && (f[0].kind() == op::And || f[0].kind() == op::Or || f[0].kind() == op::G || f[0].kind() == op::F)) {
	// 	node->value = f[0].kind();
	// 	for (formula not1 : f[0]) 
	// 	{
	// 		formula Not = formula::Not(not1);
	// 		TreeNode* child = createTree(Not);
	// 		node->children.push_back(child);
	// 	}
	// }
	else if (node->value == op::G && f[0].kind() == op::And)
    {
		node->value = op::And;
		for (formula conjunct : f[0])
		{
			formula gConjunct = formula::G(conjunct);
			TreeNode* child = createTree(gConjunct);
			node->children.push_back(child);
		}
    } else if (node->value == op::F && f[0].kind() == op::Or)
    {
		node->value = op::Or;
		for (formula sum : f[0])
		{
			formula fSum = formula::F(sum);
			TreeNode* child = createTree(fSum);
			node->children.push_back(child);
		}
    } else if (node -> value == op::Implies) {
		node->value = op::Or;
		formula c1 = formula::Not(f[0]);
		TreeNode* child1 = createTree(c1);
		formula c2 = f[1];
		TreeNode* child2 = createTree(c2);
		node->children.push_back(child1);
		node->children.push_back(child2);
	}
	
	else
    {
		twa_graph_ptr aut = minimize_explicit(trans_formula(f, dict, opt->_num_ap_for_mona));
        //cout << aut->num_states() << endl;
        dfwa_pair* pair = new dfwa_pair(aut, aut->num_states(), true, f);

		node->dfa = pair;
        //node->leafValue = f;
		sizes++;
    }

    return node;
}

void deleteTree(TreeNode* node)
{
    if (node == nullptr)
        return;

    for (auto child : node->children)
        deleteTree(child);

    delete node;
}

void combineNot(TreeNode* node) {
	//node->leafValue = formula::Not((node->children).front()->leafValue);
	//twa_graph_ptr childDFACompliment = spot::complement(((node->children).front()->dfa)->_twa);
	//node->dfa = new dfwa_pair(childDFACompliment, childDFACompliment->num_states(), true, node->leafValue);
}

void combineAnd(TreeNode* node) {
	priority_queue<dfwa_pair, std::vector<dfwa_pair>, GreaterThanByDfwaSize> childrenSorted;
	for (TreeNode* child : node->children) {
		childrenSorted.push(*(child->dfa));
		delete child;
	}
	while (childrenSorted.size() > 1) {
		dfwa_pair first = childrenSorted.top();
        childrenSorted.pop();
        dfwa_pair second = childrenSorted.top();
        childrenSorted.pop();

        formula result_formula = formula::And({first._formula, second._formula});
		twa_graph_ptr A = first._twa;
		twa_graph_ptr B = second._twa;

		twa_graph_ptr P = spot::product(A, B);
		P = minimize_explicit(P);
		dfwa_pair pair = dfwa_pair(P, P->num_states(), true, result_formula);
		pair._num_propduct = 1;
		childrenSorted.push(pair);
	}

	node->dfa = new dfwa_pair(childrenSorted.top());
	//node->leafValue = node->dfa->_formula;
}

void combineOr(TreeNode* node) {
	priority_queue<dfwa_pair, std::vector<dfwa_pair>, GreaterThanByDfwaSize> childrenSorted;
	for (TreeNode* child : node->children) {
		childrenSorted.push(*(child->dfa));
		delete child;
	}
	while (childrenSorted.size() > 1) {
		dfwa_pair first = childrenSorted.top();
        childrenSorted.pop();
        dfwa_pair second = childrenSorted.top();
        childrenSorted.pop();

        formula result_formula = formula::Or({first._formula, second._formula});
		twa_graph_ptr A = first._twa;
		twa_graph_ptr B = second._twa;

		twa_graph_ptr P = spot::product_or(A, B);
        
		P = minimize_explicit(P);
		dfwa_pair pair = dfwa_pair(P, P->num_states(), true, result_formula);
		pair._num_propduct = 1;
		childrenSorted.push(pair);
	}

	node->dfa = new dfwa_pair(childrenSorted.top());
	//node->leafValue = node->dfa->_formula;
}

void combineTree(TreeNode* node)
{
    if (node == nullptr || node->children.empty())
        return;

	for (TreeNode* child : node->children)
		combineTree(child);

	if (node->value == spot::op::Not)
	{
		combineNot(node);
	}
	else if (node->value == spot::op::And)
	{
		combineAnd(node);
	}
	else if (node->value == spot::op::Or)
	{
		combineOr(node);
	}
} 

void printTree(TreeNode* node, int indent = 0)
{
    if (node == nullptr)
        return;

    std::string indentStr(indent, ' ');

    if (node->value == op::And || node->value == op::Or 
	//node->value == op::Not
	)
    {
        std::cout << indentStr << "Operator: " << opToString(node->value) << std::endl;
    }
    else
    {
        std::cout << indentStr << "Formula: " << node->dfa->_formula << std::endl;
    }

    for (auto child : node->children)
        printTree(child, indent + 4);
}


int main(int argc, char** argv)
{
    opt_t o;
    opt = &o;
    parse_opt(argc, argv);
	
    ifstream ltlfile;
	ltlfile.open(opt->_ltlfile_name);
    string line;
    clock_t c_start = clock();
    formula input_f;
	// cout << ltlfile << endl;
	//cout << "Starting the decomposition phase" << endl;

    getline(ltlfile, line);
    //cout << "formula: " << line << endl;
    auto pf1 = spot::parse_infix_psl(line.c_str());
    if (pf1.format_errors(std::cerr))
    {
        std::cerr << "Error: " << line << std::endl;
        return -1;
    }
    // formula 
    input_f = pf1.f;

    TreeNode* tree = createTree(input_f);
	//printTree(tree);
	clock_t c_end = clock();
	cout << "Decomp Runtime: " << 1000.0 * (c_end - c_start)/CLOCKS_PER_SEC << "ms ..." << endl;
	cout << "Breakdown: " << sizes << endl;
    ltlfile.close();

    //cout << "Starting the composition phase" << endl;

    combineTree(tree);

	c_end = clock();
	//cout << tree->leafValue << endl;

	dfwa_pair* pair = tree->dfa;

    // cout << "Finished constructing minimal dfa in "
    // 	<< 1000.0 * (c_end - c_start)/CLOCKS_PER_SEC << "ms ..." << endl;
	// cout << "Number of states (or nodes) is: " << pair->_num_states << endl;
	pair->_twa = minimize_explicit(pair->_twa);
	// cout << "Final result (or number of nodes): " << pair->_twa->num_states() << endl;

	cout << "Breakdown: " << sizes << endl;
	cout << "Runtime: " << 1000.0 * (c_end - c_start)/CLOCKS_PER_SEC << "ms ..." << endl;
	cout << "Min States: " << pair->_twa->num_states() << endl;


}