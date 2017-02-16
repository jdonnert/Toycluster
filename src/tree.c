#include "globals.h"

#define NODES_PER_PARTICLE 2.5

struct Tree_Node {
	uint32_t Bitfield; 	// bit 0-5:level, 6-8:key, 9:local, 10:top, 11-31:free
	int DNext;		   	// Distance to the next node; or particle -DNext-1
	float Pos[3];		// Node Center
	int Npart;			// Number of particles in node
	float Size;
} *Tree;

int NNodes = 0;
int Max_Nodes = 0;

static inline int key_fragment(const int node);
static inline void add_particle_to_node(const int ipart, const int node);
static inline bool particle_is_inside_node(const peanoKey key, const int lvl,		const int node);
static inline void create_node_from_particle(const int ipart,const int parent,
		const peanoKey key, const int lvl, int *NNodes);
void Tree_Build();
void gravity_tree_init();
int Level(const int node); 

int Find_ngb_tree(const int ipart, const float hsml, int ngblist[NGBMAX])
{    
	const float boxsize =  Param.Boxsize;
    const float boxhalf = Param.Boxsize * 0.5;
	const float pos_i[3] = {P[ipart].Pos[0],P[ipart].Pos[1],P[ipart].Pos[2]};

	int node = 1;

	int ngbcnt = 0;

	for (;;) {
    	
		float dx = fabs(pos_i[0] - Tree[node].Pos[0]);
        float dy = fabs(pos_i[1] - Tree[node].Pos[1]);
        float dz = fabs(pos_i[2] - Tree[node].Pos[2]);
			
        if (dx > boxhalf)
	 	  	dx -= boxsize;
		
        if (dy > boxhalf)
	 	  	dy -= boxsize;
		
        if (dz > boxhalf)
	 	  	dz -= boxsize;

		float dl = 0.5 * sqrt3 * Tree[node].Size + hsml;
		
		if (dx*dx + dy*dy + dz*dz < dl*dl) {

			if (Tree[node].DNext < 0) {
			
				int first = -(Tree[node].DNext + 1);
				int last = first + Tree[node].Npart;

				for (int jpart = first; jpart < last; jpart++) { 
				
					float dx = fabs(pos_i[0] - P[jpart].Pos[0]);
           			float dy = fabs(pos_i[1] - P[jpart].Pos[1]);
           			float dz = fabs(pos_i[2] - P[jpart].Pos[2]);

        			if (dx > boxhalf)
				 	  	dx -= boxsize;

	        		if (dy > boxhalf)
				 	  	dy -= boxsize;

			        if (dz > boxhalf)
	 	  				dz -= boxsize;

					float dl = 0.5 * sqrt3 * Tree[node].Size + hsml;

					if (dx*dx + dy*dy + dz*dz < hsml*hsml)
						ngblist[ngbcnt++] = jpart;

					if (ngbcnt == NGBMAX)
						return ngbcnt;
				}
			}

			node++; 

			if (node >= NNodes) 
				break;
			else
				continue;
		} // if dx < dl 
				
		node += max(1, Tree[node].DNext);
		
		if (node >= NNodes) 
			break;
	}

	return ngbcnt;
}

extern float Guess_hsml(const size_t ipart, const int DesNumNgb)
{
	int node = P[ipart].Tree_Parent;

    float numDens = Tree[node].Npart / p3(Tree[node].Size);
    float size = pow( fourpithird/numDens, 1./3.);
    
	return 2 * size;
}


void Build_Tree()
{
	gravity_tree_init();
	
	const double boxsize = Param.Boxsize;

	create_node_from_particle(0, 0, 0, 0, &NNodes); // root node 

	Tree[0].Pos[0] = Tree[0].Pos[1] = Tree[0].Pos[2] = boxsize/2;

	int last_parent = 0; // last parent of last particle

	double px = P[0].Pos[0]/boxsize;
	double py = P[0].Pos[1]/boxsize;
	double pz = P[0].Pos[2]/boxsize;
		
	peanoKey last_key = Reversed_Peano_Key(px, py, pz);
	
	last_key >>= 3; 

	for (int ipart = 1; ipart < Param.Npart[0]; ipart++) {

		double px = P[ipart].Pos[0]/boxsize; 	
		double py = P[ipart].Pos[1]/boxsize; 
		double pz = P[ipart].Pos[2]/boxsize; 
		
		peanoKey key = Reversed_Peano_Key(px, py, pz);

		int node = 0; // current node
		int lvl = 0; // counts current level
		int parent = 0; // parent of current node

		while (lvl < N_PEANO_TRIPLETS) {
			
			if (particle_is_inside_node(key, lvl, node)) { // open node	
				
				if (Tree[node].Npart == 1) { // refine 
	
					Tree[node].DNext = 0;		

					create_node_from_particle(ipart-1, node, last_key, lvl+1, 
							&NNodes); // son of node

					last_key >>= 3;
				}  
				
				add_particle_to_node(ipart, node); // add ipart to node

				parent = node;

				node++; // decline into node

				lvl++;
				
				key >>= 3;

			} else { // skip node
				
				if (Tree[node].DNext == 0 || node == NNodes - 1)   
					break; // reached end of branch
				
				node += fmax(1, Tree[node].DNext);
			}
		} // while (lvl < 42)
		
		if (lvl > N_PEANO_TRIPLETS-1) {	// particles closer than PH resolution
		
			P[ipart].Tree_Parent = parent;
			
			continue; 					// tree cannot be deeper
		}

		if (Tree[node].DNext == 0) 				// set DNext for internal node
			Tree[node].DNext = NNodes - node; 	// only delta
			
		create_node_from_particle(ipart, parent, key, lvl, &NNodes); // sibling
	
		last_key = key >> 3;
		last_parent = parent;

	} // for ipart

	Tree[0].DNext = 0; 

	int stack[N_PEANO_TRIPLETS + 1] = { 0 }; 
	int lowest = 0;

	for (int i = 1; i < NNodes; i++) {
		
		int lvl = Level(i);

		while (lvl <= lowest) { // set pointers

			int node = stack[lowest];
	
			if (node > 0)
				Tree[node].DNext = i - node;

			stack[lowest] = 0;

			lowest--;
		} 
		
		if (Tree[i].DNext == 0) { // add node to stack
			
			stack[lvl] = i;
			
			lowest = lvl;
		}
		
	} // for

	return ;
}

static inline bool particle_is_inside_node(const peanoKey key, const int lvl,		const int node)
{
	int part_triplet = key & 0x7;

	int node_triplet = key_fragment(node); 

	return (node_triplet == part_triplet); 
}

static inline void create_node_from_particle(const int ipart,const int parent, 
		const peanoKey key, const int lvl, int *NNodes)
{
	const int node = (*NNodes)++;

	Assert(*NNodes < Max_Nodes, "Too many tree nodes %d \n", Max_Nodes);

	Tree[node].DNext = -ipart - 1;

	int keyfragment = (key & 0x7) << 6;

	Tree[node].Bitfield = lvl | keyfragment;

	const int sign[3] = { -1 + 2 * (P[ipart].Pos[0] > Tree[parent].Pos[0]),
	 			     	  -1 + 2 * (P[ipart].Pos[1] > Tree[parent].Pos[1]),
	 			          -1 + 2 * (P[ipart].Pos[2] > Tree[parent].Pos[2]) }; 
	
	float size = Param.Boxsize / (1 << lvl);

	Tree[node].Size = size;

	Tree[node].Pos[0] = Tree[parent].Pos[0] + sign[0] * size * 0.5;
	Tree[node].Pos[1] = Tree[parent].Pos[1] + sign[1] * size * 0.5;
	Tree[node].Pos[2] = Tree[parent].Pos[2] + sign[2] * size * 0.5;

	P[ipart].Tree_Parent = parent;

	add_particle_to_node(ipart, node); 

	return ;
}

static inline void add_particle_to_node(const int ipart, const int node)
{
	Tree[node].Npart++;
	
	return ;
}

static inline int key_fragment(const int node)
{
	const uint32_t bitmask = 7UL << 6;

	return (Tree[node].Bitfield & bitmask) >> 6; // return bit 6-8
}

int Level(const int node)
{
	return Tree[node].Bitfield & 0x3FUL; // return but 0-5
}


void gravity_tree_init()
{
	Max_Nodes = Param.Npart[0] * NODES_PER_PARTICLE;
	
	size_t nBytes = Max_Nodes * sizeof(*Tree);

	if (Tree == NULL)
		Tree = malloc(nBytes);
	
	memset(Tree, 0, nBytes);

	NNodes = 0;

	return ;
}
