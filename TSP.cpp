#include <fstream>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <vector>
#include <iostream>
#include <string>
#include <cmath>



std::vector<std::vector<std::vector<int> > > enumeration(int n, const std::string& combtype);
bool isComb(const std::vector<std::vector<int> >& comb, const std::string& combtype);
std::vector<std::vector<int> > subsets(const std::vector<int>& nodes);
std::vector<std::vector<std::vector<int> > > subsetcomb(const std::vector<std::vector<int> >& subsetsofnodes, const std::string& combtype);
void calcSubset(const std::vector<int>& nodes, std::vector<std::vector<int> >& res, std::vector<int>& subset, int index);
void calcSubsetComb(const std::vector<std::vector<int> >& subsetsofnodes, std::vector<std::vector<std::vector<int> > >& res, std::vector<std::vector<int> >& subset, int index, const std::string& combtype);
int getedge(int vector1, int vector2, int n);
std::vector<int> reverseedge(int edge_num, int n);
void combinequality(std::vector<std::vector<int> >& A, const std::vector<std::vector<std::vector<int> > >& combs, int n);
void STinequalities(std::vector<std::vector<int> >& A, int n);
std::vector<double> metric_dist(const std::vector<std::vector<double> >& coordinates, int n);
void printvec(const std::vector<int>& vec);
void printvec(const std::vector<double>& vec);
void printmatrix(const std::vector<std::vector<int> >& matrix);
void createfile(const std::vector<std::vector<int> >& A, const std::vector<double>& opt);
void finderror(const std::vector<double> & solution, const std::vector<std::vector<int> > & A);

//gives the option to start on the first or second column of the matrix (lpsolver doesnt use the index 0). 0 for the first and 1 for the second.
int startofrow = 1; 

int main(){
    int n = 6;
    // combtype determines the type of combs we want to use for our solution, allowed inputs are "subtour" "comb" "simple" "chvatal" and "blossom"
    std::string combtype = "subtour";
    time_t start, middle, end;
    double time_taken;
    std::vector<std::vector<int> > A;
    if(combtype != "subtour"){
        time(&start);
        std::vector<std::vector<std::vector<int> > > combs = enumeration(n, combtype);
        time(&middle);
        time_taken = double(middle - start);
        std::cout << "Time taken for subsets for n = " << n <<  " is : " << time_taken << " sec " << std::endl;
        combinequality(A, combs, n);
    }
    else{
        time(&middle);
    }
    STinequalities(A,n);
    time(&end);
    time_taken = double(end - middle);
    std::cout << "Time taken for " << combtype << " comb and subtour inequalities for n = " << n <<  " is : " << time_taken << " sec " << std::endl;

    //finderror({0, 0.55, 0, 0.45, 0.55, 0.45, 0, 0, 0.45, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0.55, 0, 0, 0, 0.5, 1, 0, 0.5, 0.55, 0, 0.45}, A);
    
    double bignum = 100;
    std::vector<double> cost = {0, 2, 3, 2, bignum, bignum, 2, bignum, 0, bignum, bignum, bignum, 2, 2, 3, 2};
    createfile(A, cost);
    

    return 0;
}


void finderror(const std::vector<double> & solution, const std::vector<std::vector<int> > & A){
    //this function finds all comb inequalities that are broken by a given solution
    for(int i = 0; i<A.size(); i++){
        double lhs = 0;
        std::vector<double> values = {};
        for (int j = 0; j<A[i].size()-2; j++){
            lhs += A[i][j]*solution[j];
            values.push_back(A[i][j]*solution[j]);
        }
        if(lhs > A[i][A[i].size()-1]){
            std::cout << lhs << std::endl;
            printvec(values);
            printvec(A[i]);
            std::cout << i << std::endl;
        }
    }
}

void createfile(const std::vector<std::vector<int> >& A, const std::vector<double>& opt){
    //This function transfer the LP stored in the matrix A to a mps file that can be read by external solvers
    std::ofstream MyFile("LP.mps");
    MyFile << "NAME" << "          " << "TEST" << '\n';
    MyFile << "ROWS" << '\n';
    MyFile << " " << "N" << "  " << "COST" <<'\n';
    for(int row_num = 0; row_num < A.size(); row_num++){
        if(A[row_num][A[row_num].size()-2] == 1){
            MyFile << " " << "L" << "  ";
        }
        else if(A[row_num][A[row_num].size()-2] == 3){
            MyFile << " " << "E" << "  ";
        }
        else{
            std::cout << "Invalid inequality type in line " << row_num << std::endl; 
        }
        MyFile << "LIM" << row_num << '\n';
    }
    MyFile << "COLUMNS" << '\n';
    for(int column_num = startofrow; column_num < A[0].size()-2; column_num++){
        MyFile << "    " << "X" << column_num;
            if(column_num < 10){
                MyFile << "        ";
            }
            else if(column_num > 9 and column_num < 100){
                MyFile << "       ";
            }
            else if(column_num > 99 and column_num < 1000){
                MyFile << "      ";
            }
            else{
                std::cout << "File too big to output to .mps" << std::endl;
            }
            MyFile << "COST" << "                " << opt[column_num-1+startofrow] << '\n';
        for(int row_num = 0; row_num < A.size(); row_num++){
            if(A[row_num][column_num] != 0){
                MyFile << "    " << "X" << column_num;
                if(column_num < 10){
                    MyFile << "        ";
                }
                else if(column_num > 9 and column_num < 100){
                    MyFile << "       ";
                }
                else if(column_num > 99 and column_num < 1000){
                    MyFile << "      ";
                }
                else{
                    std::cout << "File too big to output to .mps" << std::endl;
                }
                MyFile << "LIM" << row_num;
                if(row_num < 10){
                    MyFile << "                ";
                }
                else if(row_num > 9 and row_num < 100){
                    MyFile << "               ";
                }
                else if(row_num > 99 and row_num < 1000){
                    MyFile << "              ";
                }
                else if(row_num > 999 and row_num < 10000){
                    MyFile << "             ";
                }
                else if(row_num > 9999 and row_num < 100000){
                    MyFile << "            ";
                }
                else if(row_num > 99999 and row_num < 1000000){
                    MyFile << "           ";
                }
                else if(row_num > 999999 and row_num < 10000000){
                    MyFile << "          ";
                }
                else if(row_num > 9999999 and row_num < 100000000){
                    MyFile << "         ";
                }
                else{
                    std::cout << "File too big to output to .mps" << std::endl;
                }
                MyFile << A[row_num][column_num] << '\n';
            }
        }
    }
    MyFile << "RHS" <<'\n';
    for(int row_num = 0; row_num<A.size(); row_num++){
        MyFile << "    " << "RHS1" << "      " << "LIM" << row_num;
        if(row_num < 10){
            MyFile << "                ";
        }
        else if(row_num > 9 and row_num < 100){
            MyFile << "               ";
        }
        else if(row_num > 99 and row_num < 1000){
            MyFile << "              ";
        }
        else if(row_num > 999 and row_num < 10000){
            MyFile << "             ";
        }
        else if(row_num > 9999 and row_num < 100000){
            MyFile << "            ";
        }
        else if(row_num > 99999 and row_num < 1000000){
            MyFile << "           ";
        }
        else if(row_num > 999999 and row_num < 10000000){
            MyFile << "          ";
        }
        else if(row_num > 9999999 and row_num < 100000000){
            MyFile << "         ";
        }
        else{
            std::cout << "File too big to output to .mps" << std::endl;
        }
        MyFile << A[row_num][A[row_num].size()-1] << '\n';
    }

    MyFile << "BOUNDS" <<'\n';
    for(int column_num = startofrow; column_num<A[0].size()-2; column_num++){
        MyFile << " " << "LO" << " " << "BND1" << "      " << "X" << column_num;
        if(column_num < 10){
            MyFile << "        ";
        }
        else if(column_num > 9 and column_num < 100){
            MyFile << "       ";
        }
        else if(column_num > 99 and column_num < 1000){
            MyFile << "      ";
        }
        else{
            std::cout << "File too big to output to .mps" << std::endl;
        }
        MyFile << 0 << '\n';
    }

    MyFile << "ENDATA";

    MyFile.close();
}


std::vector<std::vector<std::vector<int> > > enumeration(int n, const std::string& combtype){
    //this function enumerates all comb-like structures (combs with a any number of teeth)
    std::vector<int> nodes;
    for(int node_num = 0; node_num<n; node_num++){
        nodes.push_back(node_num);
    }
    std::vector<std::vector<int> > res = subsets(nodes);
    std::vector<std::vector<std::vector<int> > > combs;
    if(combtype == "comb" or combtype == "simple" or combtype == "chvatal" or combtype == "blossom"){
        combs  = subsetcomb(res, combtype);
    }
    else{
        std::cout << "Error with type of comb" << std::endl;
    }

    return combs;
}


bool isComb(const std::vector<std::vector<int> >& comb, const std::string& combtype){
    //this function checks if a given set of subsets is a comb, excluding the condition for the number of teeth
    //this for loop checks for intersections with the head
    for(int comb_element_num = 1; comb_element_num<comb.size(); comb_element_num++){
        int sizeofcut = 0;
        bool validtoothin = false;
        bool validtoothout = false;
        if(comb[comb_element_num].size() != 2 and combtype == "blossom"){
            return false;
        }
        for(int tooth_node_num = 0;  tooth_node_num < comb[comb_element_num].size(); tooth_node_num++){
            bool validtoothouttemp = true;
            for(int head_node_num = 0; head_node_num < comb[0].size(); head_node_num++){
                if(comb[0][head_node_num] == comb[comb_element_num][tooth_node_num]){
                    validtoothin = true;
                    validtoothouttemp = false;
                    sizeofcut += 1;
                }
            }
            if(validtoothouttemp){
                validtoothout = true;
            }
        }
        if(sizeofcut != 1 and sizeofcut != comb[comb_element_num].size()-1 and combtype == "simple"){
            return false;
        }
        if(sizeofcut != 1 and combtype == "chvatal"){
            return false;
        }
        if(sizeofcut != 1 and combtype == "blossom"){
            return false;
        }
        if(validtoothin == false){
            return false;
        }
        if(validtoothout == false){
            return false;
        }
    }
    //this for loop checks for intersections between the teeth
    for(int tooth1_num = 1; tooth1_num<comb.size(); tooth1_num++){
        for(int tooth2_num = tooth1_num+1; tooth2_num<comb.size(); tooth2_num++){
            for(int tooth1_node_num = 0; tooth1_node_num < comb[tooth1_num].size(); tooth1_node_num++){
                for(int tooth2_node_num = 0; tooth2_node_num < comb[tooth2_num].size(); tooth2_node_num++){
                if(comb[tooth1_num][tooth1_node_num] == comb[tooth2_num][tooth2_node_num]){
                    return false;
                }
                }
            }
        }
    }
    return true;
}



std::vector<std::vector<int> > subsets(const std::vector<int>& nodes){
    //this is a wrapper function for calcSubsetComb
    std::vector<int> subset;
    std::vector<std::vector<int> > res;
    int index = 0;
    calcSubset(nodes, res, subset, index);
    return res;
}

std::vector<std::vector<std::vector<int> > > subsetcomb(const std::vector<std::vector<int> >& subsetsofnodes, const std::string& combtype){
    //this is a wrapper function for calcSubsetComb
    std::vector<std::vector<int> > subset;
    std::vector<std::vector<std::vector<int> > > res;
    int index = 0;
    calcSubsetComb(subsetsofnodes, res, subset, index, combtype);
    return res;
}
 
void calcSubset(const std::vector<int>& nodes, std::vector<std::vector<int> >& res, std::vector<int>& subset, int index){
    //this function adds to res a set of subsets of std::vector A
    // Add the current subset to the result list
    res.push_back(subset);
 
    // Generate subsets by recursively including and
    // excluding elements
    for (int i = index; i < nodes.size(); i++) {
        // Include the current element in the subset
        subset.push_back(nodes[i]);
 
        // Recursively generate subsets with the current
        // element included
        calcSubset(nodes, res, subset, i + 1);
 
        // Exclude the current element from the subset
        // (backtracking)
        subset.pop_back();
    }
}

void calcSubsetComb(const std::vector<std::vector<int> >& subsetsofnodes, std::vector<std::vector<std::vector<int> > >& res, std::vector<std::vector<int> >& subset, int index, const std::string& combtype){
    //this function adds to res all comb-like elements that are a subset of A (which will be itself a set of subsets)
    // Add the current subset to the result list
    res.push_back(subset);
    // Generate subsets by recursively including and
    // excluding elements
    for (int i = index; i < subsetsofnodes.size(); i++) {
            // Include the current element in the subset
            subset.push_back(subsetsofnodes[i]);
        if(isComb(subset, combtype)){
            // Recursively generate subsets with the current
            // element included
            if(subset.size() == 1){
                calcSubsetComb(subsetsofnodes, res, subset, 0, combtype);
            }
            else{
                calcSubsetComb(subsetsofnodes, res, subset, i + 1, combtype);
            }
        }
            // Exclude the current element from the subset
            // (backtracking)
            subset.pop_back();
        
    }
}


int getedge(int vector1, int vector2, int n){       
    //this function gives an enumeration for the edges which we use to assign each edge a column in the matrix
    if(vector1 > vector2){
        int temp = vector1;
        vector1 = vector2;
        vector2 = temp;
    }
    int edge = 0;
    for(int i = 0; i<vector1; i++){
        edge+= n-i-1;
    }
    edge += vector2-vector1-1+startofrow;
    return edge;
}

std::vector<int> reverseedge(int edge_num, int n){
    int node1_num = 0;
    int node2_num = edge_num;
    while(node2_num < n-node1_num-1){
        node2_num -= n-node1_num-1;
        node1_num += 1;
    }
    std::vector<int> edge;
    edge.push_back(node1_num);
    edge.push_back(node2_num);
    return edge;
}

void combinequality(std::vector<std::vector<int> >& A, const std::vector<std::vector<std::vector<int> > >& combs, int n){      
    //this function generates a corresponding combinequality for a given comb
    for(int comb_num = 0; comb_num<combs.size(); comb_num++){
        //if it is not a comb (because it has an even number of teeth or less than 3 teeth) no inequality is added
        if(combs[comb_num].size()%2 == 0 and combs[comb_num].size()>2){
            std::vector<int> ineq;
            for (int edge_num = 0; edge_num<getedge(n,n-1, n); edge_num++){
                ineq.push_back(0);
            }
            for(int comb_element_num = 0; comb_element_num<combs[comb_num].size(); comb_element_num++){
                for(int node1_num = 0; node1_num<combs[comb_num][comb_element_num].size(); node1_num++){
                    for(int node2_num = 0; node2_num<node1_num ; node2_num++){
                        ineq[getedge(combs[comb_num][comb_element_num][node2_num], combs[comb_num][comb_element_num][node1_num], n)] += 1;
                    }
                }
            }
            int b = 0;
            for(int i = 0; i<combs[comb_num].size(); i++){
                b += combs[comb_num][i].size()-1;
            }
            b += 1;
            double k = combs[comb_num].size()-1;
            b -= std::ceil(k/2);
            //the second to last entry corresponds to the <= (1) >= (2) or = (3)
            ineq.push_back(1);
            //the last entry corresponds to the right side std::vector b
            ineq.push_back(b);
            A.push_back(ineq);
        }
    }
}

void STinequalities(std::vector<std::vector<int> >& A, int n){      
    //this function generates the subtour elimination inequalities
    //this for loop generates the first set of inequality
    for(int edge_num = startofrow; edge_num<getedge(n,n-1, n); edge_num++){
        std::vector<int> ineq;
        for (int i = 0; i<getedge(n,n-1, n); i++){
            ineq.push_back(0);
        }
        ineq[edge_num] += 1;
        //the second to last entry corresponds to the <= (1) >= (2) or = (3)
        ineq.push_back(1);
        //the last entry corresponds to the right side std::vector b
        ineq.push_back(1);
        A.push_back(ineq);
    }
    //this for loop generates the second set of inequalities
    for(int node_num = 0; node_num<n; node_num++){
        std::vector<int> ineq;
        for (int i = 0; i<getedge(n,n-1, n); i++){
            ineq.push_back(0);
        }
        for(int adjacent_node_num = 0; adjacent_node_num<n; adjacent_node_num++){
            if(node_num != adjacent_node_num){
                ineq[getedge(node_num , adjacent_node_num , n)] += 1;
            }
        }
        //the second to last entry corresponds to the <= (1) >= (2) or = (3)
        ineq.push_back(3);
        //the last entry corresponds to the right side std::vector b
        ineq.push_back(2);   
        A.push_back(ineq);
    }
    std::vector<int> nodes;
    for(int node_num = 0; node_num<n; node_num++){
        nodes.push_back(node_num);
    }
    std::vector<std::vector<int> > res = subsets(nodes);
    //this for loop generates the last set of inequalities thorough enumeration, instead of using the conventional cutting plane method
    for(int subset_num = 0; subset_num<res.size(); subset_num++){    
        //we only consider subsets with more than 2 elements as the smaller ones are redundant                
        if(res[subset_num].size() > 2 and res[subset_num].size() < n){       
            std::vector<int> ineq;
            for (int i = 0; i<getedge(n,n-1, n); i++){
                ineq.push_back(0);
            }
            for(int node1_num = 0; node1_num<res[subset_num].size(); node1_num++){
                for(int node2_num = 0; node2_num<node1_num ; node2_num++){
                    ineq[getedge(res[subset_num][node2_num], res[subset_num][node1_num], n)] += 1;
                }
            }
            //the second to last entry corresponds to the <= (1) >= (2) or = (3)
            ineq.push_back(1);
            //the last entry corresponds to the right side std::vector b
            ineq.push_back(res[subset_num].size()-1);
            A.push_back(ineq);
        }
    }

}

std::vector<double> metric_dist(const std::vector<std::vector<double> >& coordinates, int n){
    std::vector<double> c;
    for(int i = 0; i<getedge(n-1,n,n); i++){
        c.push_back(0);
    }
    for(int node1_num = 0; node1_num<n; node1_num++){
        for(int node2_num = 0; node2_num<node1_num; node2_num++){
            double distance_x = coordinates[node1_num][0]-coordinates[node2_num][0];
            double distance_y = coordinates[node1_num][1]-coordinates[node2_num][1];
            double distance = sqrt(pow(distance_x,2)+pow(distance_y,2));
            c[getedge(node1_num, node2_num, n)] = distance;
        }
    }
    return c;
}


void printvec(const std::vector<int>& vec){
    for( int i = 0; i<vec.size(); i++){
        std::cout << vec[i] << " ";
    }
    std::cout << std::endl;
}

void printvec(const std::vector<double>& vec){
    for( int i = 0; i<vec.size(); i++){
        std::cout << vec[i] << " ";
    }
    std::cout << std::endl;
}

void printmatrix(const std::vector<std::vector<int> >& vec){
    for( int i = 0; i<vec.size(); i++){
        for(int j = 0; j<vec[i].size(); j++){
            std::cout << vec[i][j] << " ";
        }
        std::cout << std::endl;
    }
}