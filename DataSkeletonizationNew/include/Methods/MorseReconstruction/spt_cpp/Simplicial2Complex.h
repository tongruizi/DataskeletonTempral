#ifndef SIMPLICIAL2COMPLEX_H
#define SIMPLICIAL2COMPLEX_H


#include <unordered_set>
#include <set>
#include <stack>
#include <map>
#include <queue>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <stack>
#include <mlpack/core.hpp>
#include "Graph.h"
#include "graph_spt.h"



using namespace std;

class Simplicial2Complex;


class Simplicial2Complex
{
    // Connectivity info
    // access these by index - which works as pointer
    // Require: no duplicate edge, no duplicate triangles
    vector<Vertex> vertexList;
    vector<Edge> edgeList;
    vector<Triangle> triList;		// this could be cleared after computing Psudo-Morse function

    vector<Vertex*> sorted_vertex;
    vector<vector<int>* > e2t;
    vector<vector<int>* > v2e;

    // filtration, only for persistence computation
    vector<Simplex*> filtration;

    // stores all critical simplesx.
    unordered_set<Simplex*> criticalSet;

    // stores all gradient arrows (index)
    DiscreteVField V;

    // stores persistence pairs (index)
    PersistencePairs P;

    int MAX_DIM = 3;
    int DIM = 3;


public:
    graph_spt sp_tree;
    vector<int> ve_e_gt_sigma;
    vector<int> et_e_gt_sigma;
    vector<int> unpaired;
    vector<int> cri_e_vec;

    vector<int> output_vert;
    vector<vector<int>> output_edge;

    vector<vector<int>> output_manifold_vert;

    void build_spt(float sigma);
    void retrieve_1manifold();
    void convert_output();
    void write_output(string vertexFile, string edgeFile);
    void write_output_vtk(std::string output);
    void write_to_graph(MyGraphType & G);

    void build_spt2(float sigma);
    void retrieve_1manifold2();
    void retrieve_1manifold_simp();

    void build_spt3(float sigma);

    void bfs(int vert_name, int &BMname, float &branch_min);
    //void bfs_p(Vertex*, int &BMname, float &branch_min);
    void dfs(int vert_name, int &BMname, float &branch_min);
    void bfs_with_min_node(int vert_name, int BMname);
    void dfs_with_min_node(int vert_name, int BMname);
    void retrive_path(int vert_name, vector<int>& vPath);

    void reset_discover_state(vector<int>);
    void set_min_for_list(vector<int>, int);

    // constructor
    Simplicial2Complex();

    // modifiers
    int addVertex(Vertex & v);
    int addVertexAlternative(Vertex & v, int i);
    int addEdge(vector<int> &v);
    int addTriangle(vector<int> &v);
    void addCriticalPoint(Simplex *s);
    void removeCriticalPoint(Simplex *s);

    // access simplex
    Vertex* atV(int i)
    {
        return &vertexList[i];
    }
    Edge* atE(int i)
    {
        return &edgeList[i];
    }
    Triangle* atT(int i)
    {
        return &triList[i];
    }
    Vertex* atS(int i)
    {
        return sorted_vertex[i];
    }
    vector<Vertex>::iterator vBegin()
    {
        return vertexList.begin();
    }
    vector<Vertex>::iterator vEnd()
    {
        return vertexList.end();
    }
    vector<Edge>::iterator eBegin()
    {
        return edgeList.begin();
    }
    vector<Edge>::iterator eEnd()
    {
        return edgeList.end();
    }
    vector<Triangle>::iterator tBegin()
    {
        return triList.begin();
    }
    vector<Triangle>::iterator tEnd()
    {
        return triList.end();
    }
    unordered_set<Simplex*>::iterator cBegin()
    {
        return criticalSet.begin();
    }
    unordered_set<Simplex*>::iterator cEnd()
    {
        return criticalSet.end();
    }
    vector<Vertex*>::iterator sBegin()
    {
        return sorted_vertex.begin();
    }
    vector<Vertex*>::iterator sEnd()
    {
        return sorted_vertex.end();
    }

    // info query
    bool isCritical(Simplex *s);
    int order();

    // connectivity operations
    vector<int>* get_edge_v(int v)
    {
        return v2e[v];
    }
    vector<int>* get_triangle_e(int e)
    {
        return e2t[e];
    }
    int oppsiteVertex(int e, int t);
    bool hasEdge(int v, int e)
    {
        /*
        for(vector<int>::iterator it = v2e[v]->begin();
        	it != v2e[v]->end(); it++){
        	if ((*it) == e){
        		return true;
        	}
        }
        return false;
        */
        int* e_vert = edgeList[e].getVertices();
        if (e_vert[0] == v || e_vert[1] == v)
            return true;
        else
            return false;
    }

    int getAdjacentVertex(int v, int e)
    {
        int* vertices = edgeList[e].getVertices();
        return vertices[0] + vertices[1] - v;
    }

    int findEdge(int v1, int v2)
    {
        for(vector<int>::iterator it = v2e[v1]->begin();
                it != v2e[v1]->end(); it++)
            if (getAdjacentVertex(v1, *it) == v2)
                return *it;
    }

    // procedural functions
    void buildComplexFromFile2_BIN(string pathname);
    void buildComplexFromInfo(arma::mat & points, std::vector<double> & f, arma::Mat<int> & tri, arma::Mat<int> & edges);
    void Load_Presaved(string input, string presave);
    void updatePsuedoMorseFunction(Edge* e);
    void buildPsuedoMorseFunction();
    void buildFiltrationWithLowerStar();
    void PhatPersistence();
    void cancelPersistencePairs(double ve_delta);
    void outputArcs(string, string, double);


    // helper functions, subroutines.
    set<Simplex*>* descendingManifold(Simplex *s);			/*Requires: s is a critical simplex*/
    void flipAndTranslateVertexFunction();
    void sortVertices()
    {
        sorted_vertex.clear();
        for(int i = 0; i < vertexList.size(); ++i)
        {
            sorted_vertex.push_back(&vertexList[i]);
        }
        sort(sorted_vertex.begin(), sorted_vertex.end(), simplexPointerCompare2);
    }
    static bool simplexPointerCompare2(const Simplex *s, const Simplex *t);
    vector<Simplex*> LowerStar(int v);
    vector<Simplex*>* isCancellable(const persistencePair01&, ofstream&);
    void cancelAlongVPath(vector<Simplex*>* VPath);
    void write_presave(string presave);

    // deprecated functions
    /*
    Vertex getVertex(int position);
    Edge getEdge(int position);
    Triangle getTriangle(int position);
    void outputComplex(string pathname);
    void buildComplexFromFile(string pathname);		// deprecated
    void buildComplexFromFile2(string pathname);	// deprecated
    void buildRipsComplex(double radius, double eps);
    DiscreteVField* getDiscreteVField(){
    	return &V;
    }
    */
};

void Simplicial2Complex::build_spt(float sigma)
{
    //also build ve_e_gt_sigma, et, unpair
    ve_e_gt_sigma.clear();
//    vector<vertex* > alloc_new;
    for(auto pp = P.msBegin(); pp != P.msEnd(); ++pp)
    {
        //add the adjacent vertices of the edge
        Edge *cri_e = atE(pp->saddle);
        int *e_vert = cri_e->getVertices();
        //assume not sorted
        Vertex *v1 = atV(e_vert[0]);
        Vertex *v2 = atV(e_vert[1]);
//        int v1_ind = v1->getoriPosition();
//        int v2_ind = v2->getoriPosition();
        int v1_ind = e_vert[0];
        int v2_ind = e_vert[1];
        float v1_val = v1->getVPosition();
        float v2_val = v2->getVPosition();

//        vertex *sp_v1 = new vertex(v1_ind, v1_val);
//        vertex *sp_v2 = new vertex(v2_ind, v2_val);
        vertex sp_v1 = vertex(v1_ind, v1_val);
        vertex sp_v2 = vertex(v2_ind, v2_val);
//        alloc_new.push_back(sp_v1);
//        alloc_new.push_back(sp_v2);
        sp_tree.add_vertex(sp_v1);
        sp_tree.add_vertex(sp_v2);

        if(pp->persistence>sigma)
        {
            //add to critical edges
            ve_e_gt_sigma.push_back(pp->saddle);
        }
        else
        {
            //add to the sp_tree
            sp_tree.add_edge(v1_ind,v2_ind);
        }
    }

//    for(int i=0;i<10;i++) {
//
//        auto test_it = sp_tree.vertices.find(i);
//        vertex vv = test_it->second;
//        cout<<vv.name<<" "<<vv.func_v<<" "<<vv.comp_min_ind<<" "<<endl;
//    }

    //sp_tree.print_graph();

//    for(auto it = alloc_new.begin();it!=alloc_new.end();it++){
//        delete [] (*it);
//    }

    et_e_gt_sigma.clear();
    for(auto pp = P.smBegin(); pp != P.smEnd(); ++pp)
    {
        if(pp->persistence>sigma)
        {
            et_e_gt_sigma.push_back(pp->saddle);
        }
    }

    unpaired.clear();
    for(auto it=edgeList.begin(); it!=edgeList.end(); it++)
    {
        if((*it).critical_type==0)
        {
            unpaired.push_back((*it).getEPosition());
        }
    }

    cout<<"ve:"<<ve_e_gt_sigma.size()<<" ";
    cout<<"et:"<<et_e_gt_sigma.size()<<" ";
    cout<<"unpair:"<<unpaired.size()<<" "<<endl;

    sp_tree.print_size();

}

void print_neighbor(Vertex* v)
{
    for(auto it=v->neighbors.begin(); it!=v->neighbors.end(); it++)
    {
        cout<<*it<<" ";
    }
    cout<<'\n';
}

void Simplicial2Complex::build_spt2(float sigma)
{
    ve_e_gt_sigma.clear();
    for(auto pp = P.msBegin(); pp != P.msEnd(); ++pp)
    {
        if (pp->persistence > sigma)
        {
            //add to critical edges
            ve_e_gt_sigma.push_back(pp->saddle);
            //delete in the simplicial complex
            Edge *cri_e = atE(pp->saddle);
            int *e_vert = cri_e->getVertices();
            Vertex *v0 = atV(e_vert[0]);
            Vertex *v1 = atV(e_vert[1]);
//            print_neighbor(v0);
            v0->neighbors.erase(v0->neighbors.find(e_vert[1]));
//            cout<<"delete "<<e_vert[1]<<endl;
            v1->neighbors.erase(v1->neighbors.find(e_vert[0]));
//            print_neighbor(v0);
//            cout<<"*****************"<<endl;

        }
    }

    et_e_gt_sigma.clear();
    for(auto pp = P.smBegin(); pp != P.smEnd(); ++pp)
    {
        if(pp->persistence>sigma)
        {
            et_e_gt_sigma.push_back(pp->saddle);
        }

        Edge *cri_e = atE(pp->saddle);
        int *e_vert = cri_e->getVertices();
        Vertex *v0 = atV(e_vert[0]);
        Vertex *v1 = atV(e_vert[1]);
        v0->neighbors.erase(v0->neighbors.find(e_vert[1]));
        v1->neighbors.erase(v1->neighbors.find(e_vert[0]));

    }
}

void Simplicial2Complex::build_spt3(float sigma)
{
    ve_e_gt_sigma.clear();
    for(auto pp = P.msBegin(); pp != P.msEnd(); ++pp)
    {

        if(pp->persistence>sigma)
        {
            //add to critical edges
            ve_e_gt_sigma.push_back(pp->saddle);

        }
        else
        {
            //add edge
            Edge *cri_e = atE(pp->saddle);
            int *e_vert = cri_e->getVertices();
            Vertex *v0 = atV(e_vert[0]);
            Vertex *v1 = atV(e_vert[1]);
            v0->neighbors.insert(e_vert[1]);
            v1->neighbors.insert(e_vert[0]);
        }
    }

    et_e_gt_sigma.clear();
    for(auto pp = P.smBegin(); pp != P.smEnd(); ++pp)
    {
        if(pp->persistence>sigma)
        {
            et_e_gt_sigma.push_back(pp->saddle);
        }
    }

    unpaired.clear();

    cout<<"ve:"<<ve_e_gt_sigma.size()<<" ";
    cout<<"et:"<<et_e_gt_sigma.size()<<" ";
    cout<<"unpair:"<<unpaired.size()<<" "<<endl;
}

void Simplicial2Complex::reset_discover_state(vector<int> node_list)
{
    for (auto it = node_list.begin(); it != node_list.end(); it++)
    {
        Vertex *node = atV(*it);
        node->discovered = false;
    }
}

void Simplicial2Complex::set_min_for_list(vector<int> node_list, int min_ind)
{
    for (auto it = node_list.begin(); it != node_list.end(); it++)
    {
        Vertex *node = atV(*it);
        node->comp_min_ind = min_ind;
    }
}

void Simplicial2Complex::bfs(int vert_name, int &BMname, float &branch_min)
{
    queue<int> q;
    vector<int> visited_list;
    q.push(vert_name);
    visited_list.push_back(vert_name);

    Vertex *vert = atV(vert_name);
    //vert->distance = 0;
    vert->pre = -1;
    branch_min = vert->getVPosition();
    BMname = vert_name;

    while (!q.empty())
    {
        int u = q.front();
        q.pop();
        Vertex *node_u = atV(u);
        if(!(node_u->discovered))
        {
            node_u->discovered = true;
            visited_list.push_back(u);
        }

        if (node_u->getVPosition() <= branch_min)
        {
            branch_min = node_u->getVPosition();
            BMname = node_u->getoriPosition();
        }

        for (auto it = node_u->neighbors.begin(); it != node_u->neighbors.end(); it++)
        {
            int v = *it;
            Vertex *node_v = atV(v);
            if (!(node_v->discovered))
            {
                node_v->discovered = true;
                visited_list.push_back(v);
                //node_v->distance = node_u->distance + 1;
                node_v->pre = node_u->getoriPosition();
                q.push(v);
            }
        }
    }
    reset_discover_state(visited_list);
    set_min_for_list(visited_list, BMname);
}

//void Simplicial2Complex::bfs_p(Vertex *vert, int &BMname, float &branch_min){
//    queue<Vertex*> q;
//    q.push(vert);
//    vector<Vertex*> visited_list;
//    int vert_name = vert->getoriPosition();
//    visited_list.push_back(vert);
//
//    //Vertex *vert = atV(vert_name);
//    //vert->distance = 0;
//    vert->pre2 = NULL;
//    branch_min = vert->getVPosition();
//    BMname = vert_name;
//
//    while (!q.empty()) {
//        Vertex *node_u = q.front();
//        int u = node_u->getoriPosition();
//        float nodeu_value = node_u->getVPosition();
//        q.pop();
//        //Vertex *node_u = atV(u);
//        if(!(node_u->discovered)) {
//            node_u->discovered = true;
//            visited_list.push_back(node_u);
//        }
//
//        if (nodeu_value <= branch_min) {
//            branch_min = nodeu_value;
//            BMname = u;
//        }
//
//        for (auto it = node_u->neighbors.begin(); it != node_u->neighbors.end(); it++) {
//            int v = *it;
//            Vertex *node_v = atV(v);
//            if (!(node_v->discovered)) {
//                node_v->discovered = true;
//                visited_list.push_back(v);
//                //node_v->distance = node_u->distance + 1;
//                node_v->pre = node_u->getoriPosition();
//                q.push(v);
//            }
//        }
//    }
//    reset_discover_state(visited_list);
//    set_min_for_list(visited_list, BMname);
//}

void Simplicial2Complex::dfs(int vert_name, int &BMname, float &branch_min)
{
    stack<int> s;
    vector<int> visited_list;
    s.push(vert_name);
    visited_list.push_back(vert_name);

    Vertex *vert = atV(vert_name);
    //vert->distance = 0;
    vert->pre = -1;
    branch_min = vert->getVPosition();
    BMname = vert_name;

    while (!s.empty())
    {
        int u = s.top();
        s.pop();
        Vertex *node_u = atV(u);
        if(!(node_u->discovered))
        {
            node_u->discovered = true;
            visited_list.push_back(u);
        }

        if (node_u->getVPosition() <= branch_min)
        {
            branch_min = node_u->getVPosition();
            BMname = node_u->getoriPosition();
        }

        for (auto it = node_u->neighbors.begin(); it != node_u->neighbors.end(); it++)
        {
            int v = *it;
            Vertex *node_v = atV(v);
            if (!(node_v->discovered))
            {
                node_v->discovered = true;
                visited_list.push_back(v);
                //node_v->distance = node_u->distance + 1;
                node_v->pre = node_u->getoriPosition();
                s.push(v);
            }
        }
    }
    reset_discover_state(visited_list);
    set_min_for_list(visited_list, BMname);
}

void Simplicial2Complex::bfs_with_min_node(int vert_name, int BMname)
{
    queue<int> q;
    vector<int> visited_list;
    q.push(vert_name);
    visited_list.push_back(vert_name);

    Vertex *vert = atV(vert_name);
    //vert->distance = 0;
    vert->pre = -1;

    while (!q.empty())
    {
        int u = q.front();
        q.pop();
        Vertex *node_u = atV(u);
        if(!(node_u->discovered))
        {
            node_u->discovered = true;
            visited_list.push_back(u);
        }

//        if (node_u->func_v <= branch_min) {
//            branch_min = node_u->func_v;
//            BMname = node_u->name;
//        }

        //Iterate over node_u neighbor
        for (auto it = node_u->neighbors.begin(); it != node_u->neighbors.end(); it++)
        {
            int v = *it;
            Vertex *node_v = atV(v);
            if (!(node_v->discovered))
            {
                node_v->discovered = true;
                visited_list.push_back(v);
                //node_v->distance = node_u->distance + 1;
                node_v->pre = node_u->getoriPosition();
                if (node_v->getoriPosition() == BMname)
                {
                    reset_discover_state(visited_list);
                    return;
                }
                q.push(v);
            }
        }
    }
    reset_discover_state(visited_list);
}

void Simplicial2Complex::dfs_with_min_node(int vert_name, int BMname)
{
    stack<int> s;
    vector<int> visited_list;
    s.push(vert_name);
    visited_list.push_back(vert_name);

    Vertex *vert = atV(vert_name);
    //vert->distance = 0;
    vert->pre = -1;

    while (!s.empty())
    {
        int u = s.top();
        s.pop();
        Vertex *node_u = atV(u);
        if(!(node_u->discovered))
        {
            node_u->discovered = true;
            visited_list.push_back(u);
        }

//        if (node_u->func_v <= branch_min) {
//            branch_min = node_u->func_v;
//            BMname = node_u->name;
//        }

        //Iterate over node_u neighbor
        for (auto it = node_u->neighbors.begin(); it != node_u->neighbors.end(); it++)
        {
            int v = *it;
            Vertex *node_v = atV(v);
            if (!(node_v->discovered))
            {
                node_v->discovered = true;
                visited_list.push_back(v);
                //node_v->distance = node_u->distance + 1;
                node_v->pre = node_u->getoriPosition();
                if (node_v->getoriPosition() == BMname)
                {
                    reset_discover_state(visited_list);
                    return;
                }
                s.push(v);
            }
        }
    }
    reset_discover_state(visited_list);
}

void Simplicial2Complex::retrive_path(int vert_name, vector<int>& vPath)
{
    vPath.clear();
    Vertex *vert = atV(vert_name);
    while(vert->pre!=-1)
    {
        vPath.push_back(vert->getoriPosition());
        int vpre_name = vert->pre;
        vert=atV(vpre_name);
    }
    vPath.push_back(vert->getoriPosition());
}

void Simplicial2Complex::retrieve_1manifold2()
{
    output_manifold_vert.clear();

    cri_e_vec.clear();
    cri_e_vec.insert(cri_e_vec.end(),ve_e_gt_sigma.begin(),ve_e_gt_sigma.end());
    cri_e_vec.insert(cri_e_vec.end(),et_e_gt_sigma.begin(),et_e_gt_sigma.end());

    cout<<"cri_e: "<<cri_e_vec.size()<<endl;
    int n_know_min=0;
    int n_not_know=0;
    for(auto it=cri_e_vec.begin(); it!=cri_e_vec.end(); it++)
    {
        int e_ind = *it;
        Edge *cri_e = atE(e_ind);
        int *e_vert = cri_e->getVertices();

        for(int i=0; i<2; i++)
        {
            int v_ind = e_vert[i];
            //BFS from v_ind
            Vertex *v = atV(v_ind);
            vector<int> vPath;
            vPath.clear();
            if(v->comp_min_ind==-1)
            {
                //not know branch min
                int minInd;
                float minV;
                //bfs(v_ind,minInd,minV);
                bfs(v_ind,minInd,minV);
                retrive_path(minInd,vPath);
                n_not_know++;
            }
            else
            {
                //know branch min
                int destInd = v->comp_min_ind;
                //bfs_with_min_node(v_ind,destInd);
                dfs_with_min_node(v_ind,destInd);
                retrive_path(destInd,vPath);
                n_know_min++;
            }
            output_manifold_vert.push_back(vPath);
        }
    }
    cout<<"#know_min:"<<n_know_min<<endl;
    cout<<"#not know min:"<<n_not_know<<endl;
//    for(auto it1=output_manifold_vert.begin();it1!=output_manifold_vert.end();it1++){
//        vector<int> path = *it1;
//        for(auto it2=path.begin();it2!=path.end();it2++){
//            cout<<*it2<<" ";
//        }
//        cout<<'\n';
//    }
}

void Simplicial2Complex::retrieve_1manifold_simp()
{
    output_manifold_vert.clear();

    cri_e_vec.clear();
    cri_e_vec.insert(cri_e_vec.end(),ve_e_gt_sigma.begin(),ve_e_gt_sigma.end());
    cri_e_vec.insert(cri_e_vec.end(),et_e_gt_sigma.begin(),et_e_gt_sigma.end());

    cout<<"cri_e: "<<cri_e_vec.size()<<endl;
    int n_know_min=0;
    int n_not_know=0;
    for(auto it=cri_e_vec.begin(); it!=cri_e_vec.end(); it++)
    {
        int e_ind = *it;
        Edge *cri_e = atE(e_ind);
        int *e_vert = cri_e->getVertices();

        for(int i=0; i<2; i++)
        {
            int v_ind = e_vert[i];
            //BFS from v_ind
            Vertex *v = atV(v_ind);
            vector<int> vPath;
            vPath.clear();
            if(v->comp_min_ind==-1)
            {
                //not know branch min
                int minInd;
                float minV;
                //bfs(v_ind,minInd,minV);
                bfs(v_ind,minInd,minV);
                bfs(minInd,minInd,minV);
                retrive_path(v_ind,vPath);
                n_not_know++;
            }
            else
            {
                //know branch min
                //int destInd = v->comp_min_ind;
                //bfs_with_min_node(v_ind,destInd);
                //dfs_with_min_node(v_ind,destInd);
                retrive_path(v_ind,vPath);
                n_know_min++;
            }
            output_manifold_vert.push_back(vPath);
        }
    }
    cout<<"#know_min:"<<n_know_min<<endl;
    cout<<"#not know min:"<<n_not_know<<endl;
//    for(auto it1=output_manifold_vert.begin();it1!=output_manifold_vert.end();it1++){
//        vector<int> path = *it1;
//        for(auto it2=path.begin();it2!=path.end();it2++){
//            cout<<*it2<<" ";
//        }
//        cout<<'\n';
//    }
}

void Simplicial2Complex::retrieve_1manifold()
{
    //put all critical edges together

    output_manifold_vert.clear();

    cri_e_vec.clear();
    cri_e_vec.insert(cri_e_vec.end(),ve_e_gt_sigma.begin(),ve_e_gt_sigma.end());
    cri_e_vec.insert(cri_e_vec.end(),et_e_gt_sigma.begin(),et_e_gt_sigma.end());
//    cri_e_vec.insert(cri_e_vec.end(),unpaired.begin(),unpaired.end());

    cout<<"cri_e: "<<cri_e_vec.size()<<endl;
    int n_know_min=0;
    int n_not_know=0;
    for(auto it=cri_e_vec.begin(); it!=cri_e_vec.end(); it++)
    {
//        Edge *cri_e = atE(pp->saddle);
//        int *e_vert = cri_e->getVertices();
//        //assume not sorted
//        Vertex *v1 = atV(e_vert[0]);
//        Vertex *v2 = atV(e_vert[1]);
////        int v1_ind = v1->getoriPosition();
////        int v2_ind = v2->getoriPosition();
//        int v1_ind = e_vert[0];
//        int v2_ind = e_vert[1];
//        float v1_val = v1->getVPosition();
//        float v2_val = v2->getVPosition();
        int e_ind = *it;
        Edge *cri_e = atE(e_ind);
        int *e_vert = cri_e->getVertices();

        //2ADJ: e_vert[0] e_vert[1]
        for(int i=0; i<2; i++)
        {
            int v_ind = e_vert[i];
            //BFS from v_ind
            vertex v = sp_tree.vertices.find(v_ind)->second;
            vector<int> vPath;
            vPath.clear();
            if(v.comp_min_ind==-1)
            {
                //not know branch min
                int minInd;
                float minV;
                sp_tree.bfs(v_ind,minInd,minV);
                sp_tree.retrive_path(minInd,vPath);
                n_not_know++;
            }
            else
            {
                //know branch min
                int destInd = v.comp_min_ind;
                sp_tree.bfs_with_min_node(v_ind,destInd);
                sp_tree.retrive_path(destInd,vPath);
                n_know_min++;
            }
            output_manifold_vert.push_back(vPath);
        }
    }
    cout<<"#know_min:"<<n_know_min<<endl;
    cout<<"#not_know_min:"<<n_not_know<<endl;

    cout<<"Retrieve done."<<endl;
//    for(auto it1=output_manifold_vert.begin();it1!=output_manifold_vert.end();it1++){
//        vector<int> path = *it1;
//        for(auto it2=path.begin();it2!=path.end();it2++){
//            cout<<*it2<<" ";
//        }
//        cout<<'\n';
//    }

}

void Simplicial2Complex::convert_output()
{
    //manifold and cri_edge
    //manifold
    output_vert.clear();
    output_edge.clear();
    map<int, int> in_out_ind;
    int vert_count=0;
    for(auto it=output_manifold_vert.begin(); it!=output_manifold_vert.end(); it++)
    {
        vector<int> one_mani = *it;
        for(int i=0; i<one_mani.size()-1; i++)
        {
            //write one edge
            int in_vind[2] = {one_mani[i],one_mani[i+1]};
            int out_vind[2] = {-1,-1};//for edge
            for(int j=0; j<2; j++)
            {
                auto map_it = in_out_ind.find(in_vind[j]);
                if(map_it!=in_out_ind.end())
                {
                    //exitst vertex
                    out_vind[j]=map_it->second;
                }
                else
                {
                    //write vertex
                    in_out_ind[in_vind[j]]=vert_count;
                    out_vind[j]=vert_count;
                    output_vert.push_back(in_vind[j]);
                    vert_count++;
                }
            }
            //write edge
            vector<int> new_e{out_vind[0],out_vind[1],-1};
            output_edge.push_back(new_e);
        }
    }
    //cri_edge ve
    for(auto it=ve_e_gt_sigma.begin(); it!=ve_e_gt_sigma.end(); it++)
    {
        int e_ind = *it;
        Edge *cri_e = atE(e_ind);
        int *e_vert = cri_e->getVertices();
        //2ADJ: e_vert[0] e_vert[1]
        int out_vind[2] = {-1,-1};//for edge
        for(int j=0; j<2; j++)
        {
            auto map_it = in_out_ind.find(e_vert[j]);
            if(map_it!=in_out_ind.end())
            {
                //exitst vertex
                out_vind[j]=map_it->second;
            }
            else
            {
                //write vertex
                in_out_ind[e_vert[j]]=vert_count;
                out_vind[j]=vert_count;
                output_vert.push_back(e_vert[j]);
                vert_count++;
            }
        }
        vector<int> new_e{out_vind[0],out_vind[1],1};
        output_edge.push_back(new_e);
    }
    //cri_edge et
    for(auto it=et_e_gt_sigma.begin(); it!=et_e_gt_sigma.end(); it++)
    {
        int e_ind = *it;
        Edge *cri_e = atE(e_ind);
        int *e_vert = cri_e->getVertices();
        //2ADJ: e_vert[0] e_vert[1]
        int out_vind[2] = {-1,-1};//for edge
        for(int j=0; j<2; j++)
        {
            auto map_it = in_out_ind.find(e_vert[j]);
            if(map_it!=in_out_ind.end())
            {
                //exitst vertex
                out_vind[j]=map_it->second;
            }
            else
            {
                //write vertex
                in_out_ind[e_vert[j]]=vert_count;
                out_vind[j]=vert_count;
                output_vert.push_back(e_vert[j]);
                vert_count++;
            }
        }
        vector<int> new_e{out_vind[0],out_vind[1],2};
        output_edge.push_back(new_e);
    }

//    //cri_edge
//    for(auto it=cri_e_vec.begin();it!=cri_e_vec.end();it++){
//        int e_ind = *it;
//        Edge *cri_e = atE(e_ind);
//        int *e_vert = cri_e->getVertices();
//        //2ADJ: e_vert[0] e_vert[1]
//        int out_vind[2] = {-1,-1};//for edge
//        for(int j=0;j<2;j++){
//            auto map_it = in_out_ind.find(e_vert[j]);
//            if(map_it!=in_out_ind.end()){
//                //exitst vertex
//                out_vind[j]=map_it->second;
//            }
//            else{
//                //write vertex
//                in_out_ind[e_vert[j]]=vert_count;
//                out_vind[j]=vert_count;
//                output_vert.push_back(e_vert[j]);
//                vert_count++;
//            }
//        }
//        vector<int> new_e{out_vind[0],out_vind[1],1};
//        output_edge.push_back(new_e);
//    }


//    for(auto it=output_edge.begin();it!=output_edge.end();it++){
//        vector<int> cur_e=*it;
//        cout<<cur_e[0]<<" "<<cur_e[1]<<endl;
//    }
    cout<<"#vert: "<<output_vert.size()<<endl;
    cout<<"#edge: "<<output_edge.size()<<endl;
}

void Simplicial2Complex::write_output(string vertexFile, string edgeFile)
{
    ofstream vFile(vertexFile);
    ofstream eFile(edgeFile);
    //write vertex
    for(auto it=output_vert.begin(); it!=output_vert.end(); it++)
    {
        int v_ind = *it;
        Vertex *v = atV(v_ind);
        for(int i=0; i<DIM; i++)
        {
            vFile<<v->getCoords()[i]<<" ";
        }
        vFile << v->getFuncValue() << " "<<endl;
    }
    vFile.close();
    //write edge
    for(auto it=output_edge.begin(); it!=output_edge.end(); it++)
    {
        vector<int> new_edge = *it;
        eFile<<new_edge[0]<<" "<<new_edge[1]<<" "<<new_edge[2]<<endl;
    }
    eFile.close();
}

void Simplicial2Complex::write_output_vtk(std::string output)
{
    std::ofstream mystream;
    mystream.open(output);
    mystream << "# vtk DataFile Version 1.0\n";
    mystream << "3D triangulation data\n";
    mystream << "ASCII\n";

    mystream << std::endl;
    mystream << "DATASET POLYDATA\n";
    mystream << "POINTS " << output_vert.size() << " float\n";
    for(auto it=output_vert.begin(); it!=output_vert.end(); it++)
    {
        int v_ind = *it;
        Vertex *v = atV(v_ind);
        for(int i=0; i<DIM; i++)
        {
            mystream<<v->getCoords()[i];
            if (i < DIM-1)
            {
                mystream << " ";
            }
        }
        mystream << std::endl;
    }
    mystream << "LINES " << (output_edge.size()) << " " << (output_edge.size())*3 << std::endl;

    for(auto it=output_edge.begin(); it!=output_edge.end(); it++)
    {
        vector<int> new_edge = *it;
        mystream<< "2 " << new_edge[0]<<" "<<new_edge[1]<< endl;
    }
//    for(auto globalit = paths.begin(); globalit != paths.end(); globalit++)
//    {
//        std::list<Point> path = *globalit;
//        for (auto pathit = path.begin(); pathit != path.end(); pathit++)
//        {
//            mystream << (*pathit) << std::endl;
//        }
//    }

}

void Simplicial2Complex::write_to_graph(MyGraphType & G)
{
    for(auto it=output_vert.begin(); it!=output_vert.end(); it++)
    {
        int v_ind = *it;
        Vertex *v = atV(v_ind);
        auto pHelper = v->getCoords();
        Point p(pHelper[0],pHelper[1],pHelper[2]);
        Graph::add_vertex(G,p);
    }
    for(auto it=output_edge.begin(); it!=output_edge.end(); it++)
    {
        vector<int> new_edge = *it;
        Graph::add_edge(G,new_edge[0],new_edge[1]);
    }




}


Simplicial2Complex::Simplicial2Complex()
{
    vertexList.clear();
    edgeList.clear();
    triList.clear();
    e2t.clear();
    v2e.clear();
    criticalSet.clear();
    sorted_vertex.clear();
    vector<vector<int>* > e2t;
    vector<vector<int>* > v2e;
    filtration.clear();
    // init V, P
}

int Simplicial2Complex::addVertexAlternative(Vertex & v, int i)
{
    v.setVposition(i);
    v.setoriposition(i);
    v.dim = 0;
    vertexList[i] = v;
    return i;

}

int Simplicial2Complex::addVertex(Vertex &  v)
{
    int position = vertexList.size();
    v.setVposition(position);
    v.setoriposition(position);
    v.dim = 0;
    vertexList.push_back(v);
    // addCriticalPoint((Simplex*) atV(position));
    return position;
}

int Simplicial2Complex::addEdge(vector<int> &v)
{
    // sorted in increasing order, use the inverse
    if (simplexPointerCompare2(atV(v[1]), atV(v[0])))
    {
        swap(v[0], v[1]);
    }
    int v1 = v[1];
    int v2 = v[0];

//    Vertex *vert1 = atV(v1);
//    Vertex *vert2 = atV(v2);
//    vert1->neighbors.insert(v2);
//    vert2->neighbors.insert(v1);

    int position = edgeList.size();

    // insert edge to v2e - done outside.

    Edge e(v1, v2);
    e.setEposition(position);
    e.critical_type = 0;
    Vertex* vp[2];
    vp[0] = atV(v1);
    vp[1] = atV(v2);
    e.set_vp(vp);

    edgeList.push_back(e);
    // addCriticalPoint((Simplex*) &edgeList[position]);
    return position;
}

int Simplicial2Complex::addTriangle(vector<int> &v)
{
    if (simplexPointerCompare2(atV(v[1]), atV(v[0])))
    {
        swap(v[1], v[0]);
    }
    if (simplexPointerCompare2(atV(v[2]), atV(v[0])))
    {
        swap(v[2], v[0]);
    }
    if (simplexPointerCompare2(atV(v[2]), atV(v[1])))
    {
        swap(v[2], v[1]);
    }
    int v1 = v[2];
    int v2 = v[1];
    int v3 = v[0];
    int position = triList.size();

    int e1 = findEdge(v1, v2);
    int e2 = findEdge(v1, v3);
    int e3 = findEdge(v2, v3);

    vector<int> vlist({v1, v2, v3});
    vector<int> elist({e1, e2, e3});
    Triangle t(vlist, elist);

    // addtriangle to e2t
    e2t[e1]->push_back(position);
    e2t[e2]->push_back(position);
    e2t[e3]->push_back(position);

    t.setTposition(position);
    Vertex* vp[3];
    vp[0] = atV(v1);
    vp[1] = atV(v2);
    vp[2] = atV(v3);
    t.set_vp(vp);

    this->triList.push_back(t);
    // addCriticalPoint((Simplex*) &triList[position]);
    return position;
}

void Simplicial2Complex::addCriticalPoint(Simplex *s)
{
    criticalSet.insert(s);
}

void Simplicial2Complex::removeCriticalPoint(Simplex *s)
{
    criticalSet.erase(s);
}

bool Simplicial2Complex::isCritical(Simplex *s)
{
    return (criticalSet.count(s) > 0);
}

int Simplicial2Complex::order()
{
    return vertexList.size() + edgeList.size() + triList.size();
}

void Simplicial2Complex::outputArcs(string vertexFile, string edgeFile, double et_delta)
{
    ofstream vFile(vertexFile);
    ofstream eFile(edgeFile);
    ofstream cri0_edge("unpair_edge.txt");

    set<Simplex*> manifolds;
    cout<< "Writing 1-stable manifold\n";

    /*
    ofstream test_o("output_simplex.txt", ios_base::trunc | ios_base::out);
    for(auto simp = cBegin(); simp != cEnd(); simp++){
    	test_o << (*simp)->dim << endl;
    }
    test_o.close();
    */
    // need reverse the function value for vertices again.
    flipAndTranslateVertexFunction();
    cout << "flipped to original vertex function values" << endl;

    int counter = 0;
    int count_cri0 = 0;
    for(unordered_set<Simplex*>::iterator it = cBegin(); it != cEnd(); it++)
    {
        Simplex *s = *it;

        if(s->dim == 1)
        {
            Edge *e = (Edge*)s;
            // For an e-t pair, if persistence is low, skip it.

            if (e->critical_type == 0)
            {
                cout<<e->critical_type<<" "<<e->persistence<<endl;
                count_cri0++;
                // int *e_vert = e->getVertices();
                // Vertex* v1 = &vertexList[e_vert[0]];
                // Vertex* v2 = &vertexList[e_vert[1]];
                // cri0_edge<<v1->getoriPosition()<<" "<<v2->getoriPosition()<<" "<<e->critical_type<<endl;
                cri0_edge<<e->getEPosition()<<endl;
            }

            if (e->critical_type == 2 && e->persistence < et_delta + EPS_compare)
                continue;
            if (e->critical_type == 1 && e->persistence < et_delta + EPS_compare)
                continue;
            // if (e->persistence < et_delta + EPS_compare) continue;

            updatePsuedoMorseFunction(e);
            double support_f = e->funcValue;
            set<Simplex*> *manifold = descendingManifold((Simplex*)e);

            for (set<Simplex*>::iterator it2 = manifold->begin(); it2 != manifold->end(); it2++)
            {
                if ((*it2)->dim == 1)
                {
                    Edge* te = (Edge*)(*it2);
                    if (te->getEval() < support_f)
                        te->setEval(support_f);
                }
                manifolds.insert(*it2);
            }
            delete manifold;
            counter++;
        }
    }
    cout<<"n critical type=0: "<<count_cri0<<endl;

    cout << "Written " << counter << "arcs\n";

    vector<Vertex*> vertices;
    vector<Edge*> edges;
    for(set<Simplex*>::iterator it = manifolds.begin(); it != manifolds.end(); it++)
    {
        Simplex* s = *it;
        if(s->dim == 0)
        {
            vertices.push_back((Vertex*)s);
        }
        else
        {
            edges.push_back((Edge*)s);
        }
    }

    // give vertices a new index - > starting from 1
    std::map<Vertex*, int> map;
    for(int i = 0; i < vertices.size(); i++)
    {
        Vertex *v = vertices[i];
        map.insert( std::pair<Vertex*,int>(v, i + 1));
        for(int j = 0; j < DIM; j++)
        {
            vFile << v->getCoords()[j] << " ";
        }
        vFile << v->getFuncValue() << " ";
        if (this->criticalSet.count((Simplex*)v) > 0)
        {
            vFile << "0";
        }
        else
        {
            vFile << "-1";
        }
        vFile << endl;
    }

    for(int i = 0; i < edges.size(); i++)
    {
        Edge *e = edges[i];
        int *e_vert = e->getVertices();
        Vertex* v1 = &vertexList[e_vert[0]];
        Vertex* v2 = &vertexList[e_vert[1]];
        eFile << map.find(v1)->second << " " << map.find(v2)->second << " ";
        if(this->criticalSet.count((Simplex*)e) > 0)
        {
            eFile << "1 ";
        }
        else
        {
            eFile << "-1 ";
        }
        eFile << e->getEval()<<" ";
        eFile << e->persistence;
        eFile << endl;
    }
}

void Simplicial2Complex::buildComplexFromFile2_BIN(string pathname)
{
    // Input filename
    ifstream file(pathname, ios::binary);

    char* int_buffer = new char[sizeof(int)];
    int* int_reader = (int*) int_buffer;
    char* double_buffer = new char[sizeof(double)];
    double* double_reader = (double*) double_buffer;

    // Read vertices.
    int numOfVertices;
    file.read(int_buffer, sizeof(int));
    numOfVertices = *int_reader;
    cout << "\tReading " << numOfVertices << "vertices" << endl;
    vertexList.reserve(numOfVertices);
    for (int i = 0; i < numOfVertices; i++)
    {
        double coords[3];
        double funcValue;
        for (int j = 0; j < 3; j++)
        {
            file.read(double_buffer, sizeof(double));
            coords[j] = *double_reader;
        }
        file.read(double_buffer, sizeof(double));
        funcValue = *double_reader;
        // funcValue = (int)(funcValue*1e5)/1.0e5;

        Vertex v(coords, funcValue);
        addVertex(v);		// all related processing moved here.
    }
    for (int i = 0; i < numOfVertices; i++)
    {
        addCriticalPoint((Simplex*) atV(i));
    }

    // Use flipped function --- maxma -> minima
    // So we can look at vertex-edge pair
    // function value is flipped back before final output.
    flipAndTranslateVertexFunction();
    cout << "\tSorting " << numOfVertices << "vertices" << endl;
    sortVertices();
    int counter = 0;
    for (auto i = sBegin(); i < sEnd(); ++i)
    {
        (*i)->setVposition(counter);
        counter++;
    }
    cout << "\tDone." << endl;

    cout << "\tPreparing adjacency graph for vertices" << endl;
    for (int i = 0; i < numOfVertices; i++)
    {
        vector<int>* adj_v = new vector<int>;
        adj_v->reserve(20);
        v2e.push_back(adj_v);
    }
    cout << "\tDone" << endl;

    // Read edges.
    int numOfEdges;
    file.read(int_buffer, sizeof(int));
    numOfEdges = *int_reader;
    cout << "\tReading " << numOfEdges << "edges" << endl;
    edgeList.reserve(numOfEdges);
    for (int i = 0; i < numOfEdges; i++)
    {
        int vIndex1, vIndex2;
        file.read(int_buffer, sizeof(int));
        vIndex1 = *int_reader;
        file.read(int_buffer, sizeof(int));
        vIndex2 = *int_reader;

        vector<int> e_vert;
        e_vert.clear();
        e_vert.push_back(vIndex1);
        e_vert.push_back(vIndex2);

        int e = addEdge(e_vert);
        // insert edge to v2e
        v2e[vIndex1]->push_back(e);
        v2e[vIndex2]->push_back(e);
    }
    for (int i = 0; i < numOfEdges; i++)
    {
        addCriticalPoint((Simplex*) atE(i));
    }
    cout << "\tDone." << endl;

    cout << "\tPreparing adjacency graph for edges" << endl;
    for (int i = 0; i < numOfEdges; i++)
    {
        vector<int>* adj_e = new vector<int>;
        adj_e->reserve(8);
        e2t.push_back(adj_e);
    }
    cout << "\tDone" << endl;

//    for (int i = 0; i < numOfVertices; i++) {
//        Vertex *v_test = atV(i);
//        for (auto it = v_test->neighbors.begin(); it != v_test->neighbors.end(); it++) {
//            cout << *it << " ";
//        }
//        cout << '\n';
//    }

    // Read triangles.
    int numOfTris;
    file.read(int_buffer, sizeof(int));
    numOfTris = *int_reader;
    cout << "\tReading " << numOfTris << "triangles" << endl;
    triList.reserve(numOfTris);
    for (int i = 0; i < numOfTris; i++)
    {
        int vIndex1, vIndex2, vIndex3;
        file.read(int_buffer, sizeof(int));
        vIndex1 = *int_reader;
        file.read(int_buffer, sizeof(int));
        vIndex2 = *int_reader;
        file.read(int_buffer, sizeof(int));
        vIndex3 = *int_reader;

        vector<int> t_vert;
        t_vert.clear();
        t_vert.push_back(vIndex1);
        t_vert.push_back(vIndex2);
        t_vert.push_back(vIndex3);

        int t = addTriangle(t_vert);
        // triangles in e2t are inserted.

    }
    for (int i = 0; i < numOfTris; i++)
    {
        addCriticalPoint((Simplex*) atT(i));
    }
    // at this point, edges triangles ues index in vertexList.
    delete int_buffer;
    delete double_buffer;

    // Debug output stream - output all simplex information in ASCII

    file.close();
    cout << "\tDone." << endl;
}

void Simplicial2Complex::buildComplexFromInfo(arma::mat & points, std::vector<double> & f, arma::Mat<int> & tri, arma::Mat<int> & edges)
{
//! Vertices:
    int numOfVertices = points.n_cols;
    vertexList.reserve(numOfVertices);

    for (int i = 0; i < points.n_cols; i++)
    {


        double coords[3];
        double functionvalue = f[i];
        for (int j = 0 ; j < 3; j++)
        {
            coords[j] = points(j,i);
        }
        Vertex v(coords,functionvalue);
        // std::cout << "Vertex Defined " << std::end;
        //addVertexAlternative(v,i);
        addVertex(v);

    }
    //! No clue what this does
    for (int i = 0; i < numOfVertices; i++)
    {
        addCriticalPoint((Simplex*) atV(i));
    }
    flipAndTranslateVertexFunction();
    cout << "\tSorting " << numOfVertices << "vertices" << endl;
    sortVertices();
    int counter = 0;
    for (auto i = sBegin(); i < sEnd(); ++i)
    {
        (*i)->setVposition(counter);
        counter++;
    }
    cout << "\tDone." << endl;
    cout << "\tPreparing adjacency graph for vertices" << endl;
    for (int i = 0; i < numOfVertices; i++)
    {
        vector<int>* adj_v = new vector<int>;
        adj_v->reserve(20);
        v2e.push_back(adj_v);
    }
    cout << "\tDone" << endl;

    //! Moving to edges:

    int numOfEdges = edges.n_cols;
    edgeList.reserve(numOfEdges);
    for (int i = 0; i < numOfEdges; i++)
    {
        int vIndex1, vIndex2;
        vIndex1 = edges(0,i);
        vIndex2 = edges(1,i);
        vector<int> e_vert;
        e_vert.clear();
        e_vert.push_back(vIndex1);
        e_vert.push_back(vIndex2);
        int e = addEdge(e_vert);
        // insert edge to v2e
        v2e[vIndex1]->push_back(e);
        v2e[vIndex2]->push_back(e);
    }
    for (int i = 0; i < numOfEdges; i++)
    {
        addCriticalPoint((Simplex*) atE(i));
    }
    cout << "\tDone." << endl;

    cout << "\tPreparing adjacency graph for edges" << endl;
    for (int i = 0; i < numOfEdges; i++)
    {
        vector<int>* adj_e = new vector<int>;
        adj_e->reserve(8);
        e2t.push_back(adj_e);
    }
    cout << "\tDone" << endl;

    //! Moving to Triangles:
    int numOfTris;
    numOfTris = tri.n_cols;
    cout << "\tReading " << numOfTris << "triangles" << endl;
    triList.reserve(numOfTris);
    for (int i = 0; i < numOfTris; i++)
    {
        int vIndex1, vIndex2, vIndex3;
        vIndex1 = tri(0,i);
        vIndex2 = tri(1,i);
        vIndex3 = tri(2,i);
        vector<int> t_vert;
        t_vert.clear();
        t_vert.push_back(vIndex1);
        t_vert.push_back(vIndex2);
        t_vert.push_back(vIndex3);

        int t = addTriangle(t_vert);
        // triangles in e2t are inserted.

    }
    for (int i = 0; i < numOfTris; i++)
    {
        addCriticalPoint((Simplex*) atT(i));
    }
    // at this point, edges triangles ues index in vertexList.
    cout << "\tDone." << endl;


}


void Simplicial2Complex::buildPsuedoMorseFunction()
{
    cout << "\t Processing "<< edgeList.size() <<" edges\n";
    for (unsigned int i = 0; i < edgeList.size(); i++)
    {
        Edge* e = atE(i);
        int* vertices = e->getVertices();
        Vertex *max = &vertexList[vertices[0]];
        Vertex *v2 = &vertexList[vertices[1]];

        if (v2->getFuncValue() > max->getFuncValue())
        {
            max = v2;
        }
        e->setFuncValue(max->getFuncValue());
    }

    cout << "\t Processing "<< triList.size() <<" triangles\n";
    for (unsigned int i = 0; i < this->triList.size(); i++)
    {
        Triangle *t = atT(i);
        int* edges = t->getEdges();
        Edge *max = atE(edges[0]);
        Edge *e2 = atE(edges[1]);
        Edge *e3 = atE(edges[2]);

        if (e2->getFuncValue() > max->getFuncValue())
        {
            max = e2;
        }
        if (e3->getFuncValue() > max->getFuncValue())
        {
            max = e3;
        }
        t->setFuncValue(max->getFuncValue());
    }
}

void Simplicial2Complex::updatePsuedoMorseFunction(Edge* e)
{
    int* vertices = e->getVertices();
    Vertex *max = &vertexList[vertices[0]];
    Vertex *v2 = &vertexList[vertices[1]];

    if (v2->getFuncValue() > max->getFuncValue())
    {
        max = v2;
    }
    e->setFuncValue(max->getFuncValue());
}

set<Simplex*>* Simplicial2Complex::descendingManifold(Simplex* s)
{
    // Discrete Vector Field: V exist here
    // this is the container for output
    set<Simplex*> *manifold;

    /*If s is a minimum, the only simplex in the descending manifold is s itself*/
    if (s->dim == 0)
    {
        manifold = new set<Simplex*>({ s });
    }
    else if (s->dim == 1)
    {
        manifold = new set<Simplex*>({ s });
        stack<Simplex*> *st = new stack<Simplex*>();
        Edge *e = (Edge*)s;
        int* e_vert = e->getVertices();

        Vertex *v1 = atV(e_vert[0]);
        Vertex *v2 = atV(e_vert[1]);

        st->push((Simplex*)v1);
        st->push((Simplex*)v2);
        manifold->insert((Simplex*)v1);
        manifold->insert((Simplex*)v2);
        while (!st->empty())
        {
            Simplex* simplex = st->top();
            st->pop();
            // simplex could either be an edge or a vertex
            if (simplex->dim == 0)
            {
                Vertex *vert = (Vertex*)simplex;
                int paired_edge = V.containsVE(vert->getoriPosition());
                // If it's a vertex, look at all incident edges and see where you can go.
                // There'll be at most one direction
                if (paired_edge >= 0)
                {
                    Edge* pe = atE(paired_edge);
                    manifold->insert((Simplex*) pe);
                    st->push((Simplex*) pe);
                }
            }
            else
            {
                // An edge can only be entered from a vertex, so that uses up its arrow
                // and it can only leave via the other vertex.

                Edge *edge = (Edge*)simplex;
                // find exit vert
                int* e_vert = edge->getVertices();
                int exitVert = e_vert[0];
                if (V.containsVE(e_vert[1]) != edge->getEPosition())
                {
                    exitVert = e_vert[1];
                }
                // insert the exit vert.
                Vertex* nvert = atV(exitVert);
                manifold->insert((Simplex*) nvert);
                st->push((Simplex*) nvert);
            }

        }

        delete st;
    }
    else
    {
        cout<<"we shouldn't be here now, SOMETHIING IS WRONG"<<endl;
    }

    return manifold;
}

void Simplicial2Complex::flipAndTranslateVertexFunction()
{
    /*Flip the function and find the maximum function value*/
    double max = 0;
    for (int i = 0; i < vertexList.size(); i++)
    {
        if (vertexList[i].getFuncValue() > max)
        {
            max = vertexList[i].getFuncValue();
        }
    }
    /*Translate by max*/
    for (int i = 0; i < vertexList.size(); i++)
    {
        double oldval = vertexList[i].getFuncValue();
        vertexList[i].setFuncValue(max - oldval);
    }
}


int Simplicial2Complex::oppsiteVertex(int e, int t)
{
    int* tVertices = triList[t].getVertices();
    int* eVertices = edgeList[e].getVertices();
    int sum = tVertices[0] + tVertices[1] + tVertices[2]
              - eVertices[0] - eVertices[1];
    return sum;
}


// tells if the 1st simplex is smaller
bool Simplicial2Complex::simplexPointerCompare2(const Simplex *s, const Simplex *t)
{
    //By function value
    double f1 = s->funcValue, f2 = t->funcValue;
    if (f1 < f2 - EPS_compare)
    {
        return true;
    }
    else if(f1 > f2 + EPS_compare)
    {
        return false;
    }
    else
    {
        //By dimension
        short int d1 = s->dim, d2 = t->dim;
        if (d1 < d2)
        {
            return true;
        }
        else if(d1 > d2)
        {
            return false;
        }
        else
        {
            //3 cases
            if (d1 == 0)
            {
                //by vPosition
                // if (((Vertex*)s)->getVPosition() == ((Vertex*)t)->getVPosition()) cout << "Caught duplicate vertex\n";
                return ((Vertex*)s)->getVPosition() < ((Vertex*)t)->getVPosition();
            }
            else if (d1 == 1)
            {
                // old version, where only vertices are compared.

                Edge* e1 = (Edge*)s;
                Edge* e2 = (Edge*)t;
                Vertex** e1v = e1->get_vp();
                Vertex** e2v = e2->get_vp();
                /*
                for(int i = 0; i < 2; ++i){
                	if (e1v[i] != e2v[i])
                		return simplexPointerCompare2((Simplex*)e1v[i], (Simplex*)e2v[i]);
                }
                cout << "Caught edge with same set of vertices\n";
                return false;
                */
                if (e1v[0] != e2v[0])
                    return simplexPointerCompare2((Simplex*)e1v[0], (Simplex*)e2v[0]);
                else
                    return e1->Grad() > e2->Grad();

            }
            else
            {
                Triangle* t1 = (Triangle*) s;
                Triangle* t2 = (Triangle*) t;

                Vertex** t1v = t1->get_vp();
                Vertex** t2v = t2->get_vp();

                //by tPosition
                for(int i = 0; i < 3; ++i)
                {
                    if (t1v[i] != t2v[i])
                        return simplexPointerCompare2((Simplex*)t1v[i], (Simplex*)t2v[i]);
                }
                cout << "Caught triangle with same set of vertices\n";
                return false;

            }
        }
    }
}


//performs cancellation
void Simplicial2Complex::cancelAlongVPath(vector<Simplex*>* VPath)
{
    // V exists
    // As long as VPath is not empty, we may assume it has at least 2 entries
    // using original vertex index.
    if (VPath->at(0)->dim == 1)
    {
        for (int i = 0; i < VPath->size(); i++)
        {
            Simplex *s = VPath->at(i);
            if (s->dim == 0)
            {
                Vertex *v = (Vertex*)s;
                if (i < VPath->size() - 1)
                {
                    Edge* e = (Edge*)(VPath->at(i + 1));
                    V.removeVE(v->getoriPosition(), e->getEPosition());
                }
                if (i > 0)
                {
                    Edge* e = (Edge*)(VPath->at(i - 1));
                    V.addVE(v->getoriPosition(), e->getEPosition());
                }
            }
        }

    }
    else if (VPath->at(0)->dim == 2)
    {
        cout << "This shouldn't happen - Vpath starting from triangle\n";
    }
    else
    {
        cout << "This shouldn't happen\n";
    }
}

/*
vector<Simplex*> Simplicial2Complex::LowerStar(int v){ // original index - OLD
	unordered_set<Simplex*> ls;
	ls.clear();

	// iterate all incident edges
	vector<int>* inci_e = get_edge_v(v);
	for (auto edge = inci_e->begin(); edge != inci_e->end(); ++edge){
		int* verts = (atE(*edge))->getVertices();
		// find the other vertex.
		int e2 = verts[0];
		if (e2 == v) e2 = verts[1];

		if (simplexPointerCompare2(atV(e2), atV(v))){
			ls.insert(atE(*edge));
			vector<int>* inci_tri = get_triangle_e(*edge);
			for (auto tri = inci_tri->begin(); tri != inci_tri->end(); ++ tri){
				int e3 = oppsiteVertex(*edge, *tri);
				if (simplexPointerCompare2(atV(e3), atV(v))){
					ls.insert(atT(*tri));
				}
			}
		}
	}

	vector<Simplex*> rtn;
	rtn.clear();
	for (auto s = ls.begin(); s != ls.end(); ++s){
		rtn.push_back(*s);
	}
	sort(rtn.begin(), rtn.end(), simplexPointerCompare2);
	return rtn;
}
*/

vector<Simplex*> Simplicial2Complex::LowerStar(int v)  // original index
{
    unordered_set<Simplex*> ls;
    ls.clear();

    // iterate all incident edges
    vector<int>* inci_e = get_edge_v(v);
    for (auto edge = inci_e->begin(); edge != inci_e->end(); ++edge)
    {
        int* verts = (atE(*edge))->getVertices();
        // find the other vertex.
        int e2 = verts[0];
        if (e2 == v)
            e2 = verts[1];

        if (simplexPointerCompare2(atV(e2), atV(v)))
        {
            ls.insert(atE(*edge));
        }
    }

    vector<Simplex*> edges;
    edges.clear();
    for (auto s = ls.begin(); s != ls.end(); ++s)
    {
        edges.push_back(*s);
    }
    sort(edges.begin(), edges.end(), simplexPointerCompare2);

    // take each sorted edge
    vector<Simplex*> rtn;
    rtn.clear();
    for(auto edge = edges.begin(); edge != edges.end(); ++edge)
    {
        Edge* e = (Edge*)(*edge);
        int* verts = e->getVertices();
        // find the other vertex.
        int e2 = verts[0];
        if (e2 == v)
            e2 = verts[1];
        int ep = e->getEPosition();
        vector<int>* inci_tri = get_triangle_e(ep);

        vector<Simplex*> triangles;
        triangles.clear();
        for (auto tri = inci_tri->begin(); tri != inci_tri->end(); ++ tri)
        {
            if (ls.count(atT(*tri)) > 0)
                continue;
            int e3 = oppsiteVertex(ep, *tri);
            if (simplexPointerCompare2(atV(e3), atV(v)))
            {
                triangles.push_back(atT(*tri));
                ls.insert(atT(*tri));
            }
        }
        sort(triangles.begin(), triangles.end(), simplexPointerCompare2);

        rtn.push_back(*edge);
        for(auto tri = triangles.begin(); tri != triangles.end(); ++tri)
        {
            rtn.push_back(*tri);
        }
    }

    return rtn;
}


void Simplicial2Complex::buildFiltrationWithLowerStar()
{
    /*Reserve space, so filtration doesn't need to resize as often*/
    filtration.reserve(order());
    unordered_set<Simplex*> dup_test;
    dup_test.clear();

    cout << "\tInserting simplicies...";
    int counter = 0;
    // use sorted vertex list
    for (vector<Vertex*>::iterator i = sBegin(); i < sEnd(); ++i)
    {
        vector<Simplex*> lower_star = LowerStar((*i)->getoriPosition());
        filtration.push_back(*i);
        (*i) ->filtrationPosition = counter ++;

        for (vector<Simplex*>::iterator j = lower_star.begin(); j < lower_star.end(); ++j)
        {
            if (dup_test.count(*j) > 0)
            {
                cout << "caught duplicate simplex";
                cout << (*j)->dim << "\n";
            }
            else
            {
                dup_test.insert(*j);
            }
            filtration.push_back(*j);
            (*j) ->filtrationPosition = counter ++;
        }
    }
    dup_test.clear();


}


void Simplicial2Complex::PhatPersistence()
{
    // generate boundary matrix
    cout << "\tInitializing boundary matrix...\n";
    cout << "\t\tMatrix size: " << this->filtration.size() << "\n";
    //phat::boundary_matrix< phat::vector_list > boundary_matrix;
    phat::boundary_matrix< phat::bit_tree_pivot_column > boundary_matrix;
    boundary_matrix.set_num_cols(this->filtration.size());

    std::vector< phat::index > temp_col;
    for (unsigned int i = 0; i < this->filtration.size(); i++)
    {
        Simplex *s = filtration[i];
        temp_col.clear();

        if (s->dim == 0)
        {
            // no boundary. no bits to set in matrix
            boundary_matrix.set_dim( i, 0 );
        }
        else if (s->dim == 1)
        {
            // edge
            boundary_matrix.set_dim( i, 1 );
            Edge *e = (Edge*)s;
            int* v = e->getVertices();
            Vertex* v1 = atV(v[0]);
            Vertex* v2 = atV(v[1]);
            if (v1->filtrationPosition < v2->filtrationPosition)
            {
                temp_col.push_back(v1->filtrationPosition);
                temp_col.push_back(v2->filtrationPosition);
            }
            else
            {
                temp_col.push_back(v2->filtrationPosition);
                temp_col.push_back(v1->filtrationPosition);
            }
        }
        else if (s->dim == 2)
        {
            // triangle
            boundary_matrix.set_dim( i, 2 );
            Triangle *t = (Triangle*)s;
            int* e = t->getEdges();
            Edge *e1 = atE(e[0]);
            Edge *e2 = atE(e[1]);
            Edge *e3 = atE(e[2]);
            Edge *edges[3] = { e1,e2,e3 };
            for(int j = 0; j < 2; j++)
            {
                for(int k = j + 1; k < 3; k++)
                {
                    if(edges[k]->filtrationPosition \
                            < edges[j]->filtrationPosition)
                    {
                        swap(edges[k], edges[j]);
                    }
                }
            }
            temp_col.push_back(edges[0]->filtrationPosition);
            temp_col.push_back(edges[1]->filtrationPosition);
            temp_col.push_back(edges[2]->filtrationPosition);
        }
        boundary_matrix.set_col( i, temp_col );
    }
    cout << "\tInitialized!\n";

    // call Phat
    cout << "\tComputing persitence pairs...\n";
    phat::persistence_pairs pairs;
    phat::compute_persistence_pairs\
    < phat::twist_reduction >( pairs, boundary_matrix );
    // pairs.sort();
    cout << "\tComputed and sorted!\n";

    // post processing: add sm-ms pairs, set critical points
    cout << "\tCounting total"<< pairs.get_num_pairs() <<" ms-pairs and sm-pairs...";
    for( phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ )
    {
        Simplex *s1 = this->filtration[pairs.get_pair( idx ).first];
        Simplex *s2 = this->filtration[pairs.get_pair( idx ).second];
        if (s1->dim == 0)
        {
            Vertex *v = (Vertex*)s1;
            Edge *e = (Edge*)s2;
            double persistence = e->funcValue - v->funcValue;
            int loc_diff = e->filtrationPosition - v->filtrationPosition;
            e->critical_type = 1;
            e->persistence = persistence;

            // use SORTED position
            persistencePair01 pp = { v->getVPosition(), e->getEPosition(), persistence, loc_diff };
            P.msinsert(pp);
            // Critical points added at the beginning
        }
        else
        {
            Edge* e = (Edge*)s1;
            Triangle* t = (Triangle*)s2;
            double persistence = t->funcValue - e->funcValue;
            int loc_diff = t->filtrationPosition - e->filtrationPosition;
            e->critical_type = 2;
            e->persistence = persistence;

            persistencePair12 pp;
            pp.saddle = e->getEPosition();
            pp.max = t->getTPosition();
            pp.persistence = persistence;
            pp.loc_diff = loc_diff;
            P.sminsert(pp);
            // Critical points added at the beginning
        }
    }
    cout << "done!\n";
}

//test cancellability for Edge - vertex pair
vector<Simplex*>* Simplicial2Complex::isCancellable(const persistencePair01& pp, ofstream& cancelData)
{
    // V exists here.
    // DiscreteVField *V = this->K->getDiscreteVField();

    Edge *e = atE(pp.saddle);
    Vertex *v = atS(pp.min);


    // if they have 0 persistence, we know they are cancellable and take this shortcut
    if (fabs(e->funcValue - v->funcValue) < EPS_compare
            && hasEdge(v->getoriPosition(),e->getEPosition()))
    {

        return new vector<Simplex*>({ (Simplex*)e, (Simplex*)v });
    }

    /*Each vector<Simplex*> is a V-path starting with Edge e and ending on a minimum Vertex
    If, after performing a DFS, there are two V-paths from e to v, then we cannot cancel them.
    If there are no V-paths from e to v, then we cannot cancel them. The only time we would be able to cancel
    is if there is exactly one V-path from e to v*/
    vector<vector<Simplex*>*> paths;

    /*An edge-vertex V-path is not capable of branching beyond the first edge, so we start with two paths*/
    vector<Simplex*> *path1 = new vector<Simplex*>({ (Simplex*)e });
    vector<Simplex*> *path2 = new vector<Simplex*>({ (Simplex*)e });
    paths.push_back(path1);
    paths.push_back(path2);

    stack<Simplex*> *st = new stack<Simplex*>();
    int * e_vert = e->getVertices();
    Vertex *v1 = atV(e_vert[0]);
    Vertex *v2 = atV(e_vert[1]);
    paths[0]->push_back((Simplex*)v1);
    paths[1]->push_back((Simplex*)v2);
    st->push((Simplex*)v1);
    st->push((Simplex*)v2);

    while (!st->empty())
    {
        Simplex* s = st->top();
        st->pop();
        if (s->dim == 0)
        {
            Vertex *vert = (Vertex*)s;
            int n_e = V.containsVE(vert->getoriPosition());
            if (n_e >= 0)
            {
                Edge* paired_edge = atE(n_e);
                /*If there is one found, add the edge to all paths ending with vert*/
                for (int i = 0; (unsigned) i < paths.size(); i++)
                {
                    vector<Simplex*> *path = paths[i];
                    if (path->back() == s)
                    {
                        path->push_back((Simplex*) paired_edge);
                    }
                }
                /*Then push the edge onto the stack*/
                st->push((Simplex*) paired_edge);
            }
            /*If there is no way out of the Vertex, then we've reached a minimum*/
        }
        else
        {
            Edge *edge = (Edge*)s;
            /*Here, we only need to check two vertices*/
            int* e_vert = edge->getVertices();

            /*We want the vertex such that the discrete gradient vector field does NOT contain
            (vertex, edge)*/
            Vertex *exitVert;

            if ((V.containsVE(e_vert[0])) != edge->getEPosition())
            {
                exitVert = atV(e_vert[0]);
            }
            else
            {
                exitVert = atV(e_vert[1]);
            }

            /*exitVert is the next on our path, so we add it to all paths ending with edge*/
            for (int i = 0; (unsigned) i < paths.size(); i++)
            {
                vector<Simplex*> *path = paths[i];
                if (path->back() == s)
                {
                    path->push_back((Simplex*)exitVert);
                }
            }
            /*Then we need to push it onto the stack*/
            st->push((Simplex*)exitVert);
        }
    }
    delete st;

    /*Now if there are none or more than one paths ending at Vertex v, the persistence
    pair is not cancellable, and we return an empty vector*/
    bool hasPath = false;
    bool has_unique = true;
    vector<Simplex*> *uniquePath = NULL;
    for (int i = 0; i < paths.size(); i++)
    {
        vector<Simplex*> *path = paths[i];
        if (path->back() == (Simplex*)v && hasPath == true)
        {
            has_unique = false;
            break;
        }
        else if (path->back() == (Simplex*)v)
        {
            hasPath = true;
            uniquePath = path;
        }
    }
    for(int i = 0; i < paths.size(); i++)
    {
        if (has_unique && paths[i]==uniquePath)
            continue;
        else
            delete paths[i];
    }
    /*If hasPath is true and the for loop didn't return as soon as it found a 2nd,
    Then uniquePath truly is unique. So we return it*/

    return uniquePath;
}

//cancels al pairs that can be cancelled
void Simplicial2Complex::cancelPersistencePairs(double ve_delta)
{
    cout << "\tSorting "<< P.mssize() << " ms-persistence pairs...\n";
    cout.flush();
    P.sortmspair();
    cout << "\tDone\n";
    cout.flush();

    ofstream persistencePairs("ve_pvalues.txt", ios_base::trunc | ios_base::out);
    // V exists here.

#if (DEBUG)
    ofstream cancelDataVE("cancelData_VE.txt", ios_base::trunc | ios_base::out);
#else
    // if not at debug mode, this file will not be created
    ofstream cancelDataVE;
#endif

    cout << "\tCancelling...\n";
    int i = 0;
    int count = 0;
    cout << "msPair: " << P.mssize() << "\tsmPair: " << P.smsize() << endl;

    for (auto pair1 = P.msBegin(); pair1 != P.msEnd(); ++pair1)
    {
        if (pair1->persistence < ve_delta + EPS_compare)
        {
            vector<Simplex*> *VPath = this->isCancellable(*pair1, cancelDataVE);

            if (VPath!=NULL)
            {
                count++;
                cancelAlongVPath(VPath);
                // I guess it won't free an Simplex*
                delete VPath;
                removeCriticalPoint(atS(pair1->min));
                removeCriticalPoint(atE(pair1->saddle));
                if (pair1->persistence > EPS_compare)
                {
                    persistencePairs << pair1->persistence << " 1\n";
                }
            }
            else
            {
                if (pair1->persistence > EPS_compare)
                {
                    persistencePairs << pair1->persistence << " -1\n";
                }
            }
        }
        else
        {
            if (pair1->persistence > EPS_compare)
            {
                persistencePairs << pair1->persistence << " -1\n";
            }
        }
        i++;

        if (i%10000==0)
        {
            cout << "\r";
            cout << "\t" << i << "/" << P.mssize() << "...";
            cout.flush();
        }
    }

    cout << "\tWriting smPair info\n";
    ofstream et_stream("et_pvalues.txt", ios_base::trunc | ios_base::out);
    P.output_sm_pair(et_stream);
    cout << "\tDone\n";

    cout << "\t-->msPair: " << count << "/" << P.mssize() <<endl;
    cout << "\tDone\n";
    cancelDataVE.close();
}


void Simplicial2Complex::Load_Presaved(string input, string presave)
{
    // almost the same as original reader, but does not sort.
    // In addition, it reads in persistence pairs.
    // Input filename
    ifstream file(input, ios::binary);

    char* int_buffer = new char[sizeof(int)];
    int* int_reader = (int*) int_buffer;
    char* double_buffer = new char[sizeof(double)];
    double* double_reader = (double*) double_buffer;

    // Read vertices.
    int numOfVertices;
    file.read(int_buffer, sizeof(int));
    numOfVertices = *int_reader;
    cout << "\tReading " << numOfVertices << "vertices" << endl;
    vertexList.reserve(numOfVertices);
    for (int i = 0; i < numOfVertices; i++)
    {
        double coords[3];
        double funcValue;
        for (int j = 0; j < 3; j++)
        {
            file.read(double_buffer, sizeof(double));
            coords[j] = *double_reader;
        }
        file.read(double_buffer, sizeof(double));
        funcValue = *double_reader;
        // funcValue = (int)(funcValue*1e5)/1.0e5;

        Vertex v(coords, funcValue);
        addVertex(v);		// all related processing moved here.
    }
    for (int i = 0; i < numOfVertices; i++)
    {
        addCriticalPoint((Simplex*) atV(i));
    }

    // Use flipped function --- maxma -> minima
    // So we can look at vertex-edge pair
    // function value is flipped back before final output.
    flipAndTranslateVertexFunction();

    // filled in sorted_vert
    // set sorted position for vertices - will be done later

    cout << "\tPreparing adjacency graph for vertices" << endl;
    for (int i = 0; i < numOfVertices; i++)
    {
        vector<int>* adj_v = new vector<int>;
        adj_v->reserve(20);
        v2e.push_back(adj_v);
    }
    cout << "\tDone" << endl;

    // Read edges.
    int numOfEdges;
    file.read(int_buffer, sizeof(int));
    numOfEdges = *int_reader;
    cout << "\tReading " << numOfEdges << "edges" << endl;
    edgeList.reserve(numOfEdges);
    for (int i = 0; i < numOfEdges; i++)
    {
        int vIndex1, vIndex2;
        file.read(int_buffer, sizeof(int));
        vIndex1 = *int_reader;
        file.read(int_buffer, sizeof(int));
        vIndex2 = *int_reader;

        vector<int> e_vert;
        e_vert.clear();
        e_vert.push_back(vIndex1);
        e_vert.push_back(vIndex2);

        int e = addEdge(e_vert);
        // insert edge to v2e
        v2e[vIndex1]->push_back(e);
        v2e[vIndex2]->push_back(e);
    }
    for (int i = 0; i < numOfEdges; i++)
    {
        addCriticalPoint((Simplex*) atE(i));
    }
    cout << "\tDone." << endl;

    cout << "\tPreparing adjacency graph for edges" << endl;
    for (int i = 0; i < numOfEdges; i++)
    {
        vector<int>* adj_e = new vector<int>;
        adj_e->reserve(8);
        e2t.push_back(adj_e);
    }
    cout << "\tDone" << endl;

    // Read triangles.
    int numOfTris;
    file.read(int_buffer, sizeof(int));
    numOfTris = *int_reader;
    cout << "\tReading " << numOfTris << "triangles" << endl;
    triList.reserve(numOfTris);
    for (int i = 0; i < numOfTris; i++)
    {
        int vIndex1, vIndex2, vIndex3;
        file.read(int_buffer, sizeof(int));
        vIndex1 = *int_reader;
        file.read(int_buffer, sizeof(int));
        vIndex2 = *int_reader;
        file.read(int_buffer, sizeof(int));
        vIndex3 = *int_reader;

        vector<int> t_vert;
        t_vert.clear();
        t_vert.push_back(vIndex1);
        t_vert.push_back(vIndex2);
        t_vert.push_back(vIndex3);

        int t = addTriangle(t_vert);
        // triangles in e2t are inserted.
    }
    for (int i = 0; i < numOfTris; i++)
    {
        addCriticalPoint((Simplex*) atT(i));
    }
    // at this point, edges triangles ues index in vertexList.
    file.close();


    // NEW part - read in Sorted Vert info
    ifstream pre_stream(presave, ios::binary);
    sorted_vertex.clear();
    sorted_vertex.reserve(numOfVertices);
    cout << "\treading " << numOfVertices << "sorted vertex indices" << endl;
    for (int i = 0; i < numOfVertices; i++)
    {
        Vertex* v = NULL;
        sorted_vertex.push_back(v);
    }
    for (int i = 0; i < numOfVertices; i++)
    {
        pre_stream.read(int_buffer, sizeof(int));
        int sorted_idx = *int_reader;
        Vertex* v = atV(i);
        v->setVposition(sorted_idx);
        sorted_vertex[sorted_idx] = atV(i);	// vposition starts with 0
    }

    // NEW - simplicial pairs
    pre_stream.read(int_buffer, sizeof(int));
    int num_ve = *int_reader;
    cout << "\treading " << num_ve << " VE pairs" << endl;
    for(int i = 0; i<num_ve; ++i)
    {
        persistencePair01 pp = PersistencePairs::read_ve_pair(pre_stream);
        P.msinsert(pp);
        // set E value
        Edge *e = atE(pp.saddle);
        e->critical_type = 1;
        e->persistence = pp.persistence;
    }

    pre_stream.read(int_buffer, sizeof(int));
    int num_et = *int_reader;
    cout << "\treading " << num_et << " ET pairs" << endl;
    for(int i = 0; i<num_et; ++i)
    {
        persistencePair12 pp = PersistencePairs::read_et_pair(pre_stream);
        P.sminsert(pp);
        // set E value
        Edge *e = atE(pp.saddle);
        e->critical_type = 2;
        e->persistence = pp.persistence;
    }

    pre_stream.close();

    delete int_buffer;
    delete double_buffer;

    // Debug output stream


    cout << "\tDone." << endl;
}

void Simplicial2Complex::write_presave(string presave)
{
    // write sorted_vertex index
    string output_name = presave + "/presave.bin";
    ofstream pre_stream(output_name, ios::binary);
    char* int_buffer = new char[sizeof(int)];
    int* int_writer = (int*) int_buffer;
    char* double_buffer = new char[sizeof(double)];
    double* double_writer = (double*) double_buffer;

    //ofstream v_order("v_order.txt", ios_base::trunc | ios_base::out);
    for (int i = 0; i < vertexList.size(); i++)
    {
        Vertex* v = atV(i);
        *int_writer = v->getVPosition();
        pre_stream.write(int_buffer, sizeof(int));

        //v_order << v->getoriPosition() << " " << v->getVPosition() <<"\n";

    }
    //v_order.close();

    // write ve pair
    int num_ve = P.mssize();
    *int_writer = num_ve;
    pre_stream.write(int_buffer, sizeof(int));

    for(auto pp = P.msBegin(); pp != P.msEnd(); ++pp)
    {
        PersistencePairs::write_ve_pair(*pp, pre_stream);
    }

//	if (DEBUG){
//		ofstream ppairs("PersistencePairs.txt", ios::binary);
//		for(auto pp = P.msBegin(); pp != P.msEnd(); ++pp){
//			PersistencePairs::write_ve_pair_debug(*pp, atV(pp->min)->getoriPosition(), ppairs);
//		}
//		ppairs.close();
//	}

    // write et pair
    int num_et = P.smsize();
    *int_writer = num_et;
    pre_stream.write(int_buffer, sizeof(int));

    for(auto pp = P.smBegin(); pp != P.smEnd(); ++pp)
    {
        PersistencePairs::write_et_pair(*pp, pre_stream);
    }


    pre_stream.close();
    cout << "Written " << vertexList.size() << "int, " << num_ve
         << "VE pair, " << num_et << "ET pair." << endl;
    delete int_buffer;
    delete double_buffer;
}

#endif // SIMPLICIAL2COMPLEX_H
