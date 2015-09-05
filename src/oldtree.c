#include "globals.h"
#include "tree.h"

#define  NODES_PER_PARTICLE 5.4

/* local prototypes */
static void init_tree();
static int idx2sign(int,int);
static int pos2idx(size_t, size_t);
static void refine(const size_t, size_t*, size_t*);

/* local vars */
struct tree_node{
    size_t down;        // To first daughter node or target part
    size_t next;        // To next subnode of father or to unkle 
    float pos[3];       // Center of node 
    int npart;          // Number of particles in node 
    float size;         // Spatial extent of node 
} *tree;

static size_t max_tree_size, treesize; // Number of nodes in tree 

extern void Build_tree()
{
    init_tree();
    
	for (size_t ipart = 0; ipart < Param.Npart[0]; ipart++) {
		
		size_t node = 0;
        size_t parent = 0;

        const float pos_i[3] = { P[ipart].Pos[0], P[ipart].Pos[1],
            P[ipart].Pos[2] };

        for (;;) {

            if (tree[node].npart > 1) { // full and not leaf, decline
                
                tree[node].npart++;
 
                size_t idx = pos2idx(ipart, node);

                parent = node;

				node = tree[node].down + idx;

            } else if (tree[node].npart == 1) { // full, but leaf, refine
            
                tree[node].npart++;
                
                refine(ipart, &parent, &node);

            } else { // empty leaf, sort in

                tree[node].down = ipart;
                
                tree[node].pos[0] = pos_i[0];
                tree[node].pos[1] = pos_i[1];
                tree[node].pos[2] = pos_i[2];

                tree[node].npart++;

                break;
            }
		} 
	}

    tree[0].npart -= 2;

	//printf("Using %zu of max. %zu treenodes for %zu particles \n" ,
    //        treesize, max_tree_size, Param.Npart[0]);

    return;
}

/* 
 * Find a first guess from the tree recursively 
 * by estimating the volume per particle from
 * the local number density and compute an hsml 
 * for NNGB neighbours 
 */

extern float Guess_hsml(const size_t ipart, const int DesNumNgb)
{
    int node = 0;
    int nextnode = tree[node].down + pos2idx(ipart,node);
	
    while (tree[nextnode].npart > DesNumNgb) {

        node = nextnode;

        nextnode = tree[node].down + pos2idx(ipart,node);
    }

    float numDens = tree[node].npart / p3(tree[node].size);
    float size = pow( fourpithird/numDens, 1./3.);

    return size;
}

/* 
 * Find neighbours closer than hsml via the tree, periodic version.
 */

extern int Find_ngb_tree(const size_t ipart, const float hsml, int *ngblist)
{
    const float boxsize =  Param.Boxsize;
    const float boxhalf = Param.Boxsize * 0.5;
    const float pos_i[3] = {P[ipart].Pos[0],P[ipart].Pos[1],P[ipart].Pos[2] };

    size_t node = 1;

    int ngbcnt = 0;


    while ((ngbcnt < NGBMAX) && node) {

        if (tree[node].npart != 0) {

            float dx = pos_i[0] - tree[node].pos[0];
            float dy = pos_i[1] - tree[node].pos[1];
            float dz = pos_i[2] - tree[node].pos[2];

            if (dx > boxhalf)
		    	dx -= boxsize;

    		if (dy > boxhalf)
	    		dy -= boxsize;

    		if (dz > boxhalf)
	    		dz -= boxsize;

			if (dx < -boxhalf) 
				dx += boxsize;

    		if (dy < -boxhalf)
	    		dy += boxsize;

			if (dz < -boxhalf)
				dz += boxsize;
            
			float dl = 2*tree[node].size + hsml;

            if (dx*dx + dy*dy + dz*dz < dl*dl) { // check for collision
                
                if (tree[node].npart > 1) { 

                    node = tree[node].down; // open
                    continue; // skip the rest
                } 
                
                ngblist[ngbcnt++] = tree[node].down; // add neighbour
            }
        } 

        node = tree[node].next; // go on
    }

  //  if (ngbcnt == NGBMAX)
   //   printf("ngblist maxed out for particle %zu \n", ipart);

	return ngbcnt;
}

static void init_tree()
{
    const int npart = Param.Npart[0];
	
    size_t nBytes =  NODES_PER_PARTICLE * npart * sizeof(struct tree_node);

    if (tree == NULL) 
        tree = Malloc(nBytes);

	tree[0].down = 1; 
	tree[0].pos[0] = tree[0].pos[1] = tree[0].pos[2] =  Param.Boxsize/2;
	tree[0].npart = 2;
	tree[0].next = 0; // the zero point is termination point
	tree[0].size = Param.Boxsize; 
    
	treesize = 9;
    max_tree_size = NODES_PER_PARTICLE * npart;

    for (int i = 0; i < 8; i++) { // do the first 8
    
		int son = 1 + i; 
    
        tree[son].npart = 0;
		tree[son].size = i;
        
        tree[son].next = son+1;

		tree[son].pos[0] = tree[son].pos[1] = tree[son].pos[2] = 0;
    }

    tree[8].next = 0;

	return;
}	


/* 
 * Refines parent node, sets size and positions from grandparent, 
 * puts old and new particle in place 
 */

static void refine(const size_t ipart, size_t *grandparent, size_t *parent)
{	
    const size_t node = *parent;
    const size_t jpart = tree[node].down;	// Remember old particle
    const float *pos_j = &tree[node].pos[0];

	tree[node].down = treesize; // Add at the end of tree 

	treesize += 8; // Add 8 more nodes 

    Assert(treesize < max_tree_size, "Tree is too large");

    int i = tree[node].size;
        
    tree[node].size = tree[*grandparent].size * 0.5; // treat parent node

    tree[node].pos[0] = tree[*grandparent].pos[0] 
			+ idx2sign(i,0) * 0.5 * tree[node].size;
    tree[node].pos[1] = tree[*grandparent].pos[1] 
			+ idx2sign(i,1) * 0.5 * tree[node].size;
    tree[node].pos[2] = tree[*grandparent].pos[2] 
			+ idx2sign(i,2) * 0.5 * tree[node].size;

	for (i = 0; i < 8; i++) { // Init new son nodes 

		int son = tree[node].down + i; 

		tree[son].npart = 0;
		tree[son].size = i;
        
        tree[son].next = son+1;

		tree[son].pos[0] = tree[son].pos[1] = tree[son].pos[2] = 0;
	}

    tree[tree[node].down + 7].next = tree[node].next; // last .next goes up
	
	int target = tree[node].down + pos2idx(jpart,node); // Treat particles
	tree[target].down = jpart;
	tree[target].npart++;

    tree[target].pos[0] = pos_j[0];
    tree[target].pos[1] = pos_j[1];
    tree[target].pos[2] = pos_j[2];

    *grandparent = node;
	*parent = tree[node].down + pos2idx(ipart,node);
	
	return ; // Return new target node for ipart
}

/* Find index of subnode from relative position */
static int pos2idx(size_t ipart,size_t node)
{
	float dx = P[ipart].Pos[0] -tree[node].pos[0];
	float dy = P[ipart].Pos[1] -tree[node].pos[1];
	float dz = P[ipart].Pos[2] -tree[node].pos[2];

	return (dx>0)  + ((dy>0) << 1) + ( (dz>0) << 2);
}

/* Convert node index to sign of component comp */
static int idx2sign(int idx, int comp)
{
	return -1+2*((idx & (1 << comp))>>comp);
}

