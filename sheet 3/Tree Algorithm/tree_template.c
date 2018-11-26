#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <time.h>


/* Let's define two types of structures, one for the particles, one for the tree nodes
 */
 
typedef struct particle_data
{
  double pos[3];
  double mass;

  double acc_tree[3];
  double acc_exact[3];
}
  particle;


typedef struct node_data
{
  double center[3];
  double len;
 
  double cm[3];
  double mass;

  struct node_data *suns[8];

  particle *p;
}
  node;


#define  MAX_POINTS               5000      /* this sets the particle number that is used */
static double opening_threshold = 0.8;      /* tree opening angle */
static double eps               = 0.001;    /* gravitational softening length */


/* Lets create, for simplicity, some static arrays which will hold our data */

#define  MAX_NODES    (5 * MAX_POINTS)

node      tree[MAX_NODES];
particle  star[MAX_POINTS];


/* this function returns a pointer to an empty tree node 
 */
node *get_empty_node(void)
{
  node *no;
  static int count_nodes = 0;

  if(count_nodes < MAX_NODES)
    {
      no = &tree[count_nodes++];
      memset(no, 0, sizeof(node));
    }
  else
    {
      printf("sorry, we are out of tree nodes.\n");
      exit(1);
    }

  return no;
}


/* this function determines the index (0-7) of the subnode in which
 * a particle falls within a given node
 */
int get_subnode_index(node *current, particle *p)
{
  int index = 0;

  if(p->pos[0] > current->center[0])
    index += 4;
  if(p->pos[1] > current->center[1])
    index += 2;
  if(p->pos[2] > current->center[2])
    index += 1;

  return index;
}


/* this function's task is to insert a new particle into a given tree node 
 */
void insert_particle(node *current, particle *pnew)
{
  node *sun;
  int i, j, k, n, p_subnode, pnew_subnode;

  if(current->p)
    {
      /* The node contains a particle already. 
	 Need to create a new set of 8 subnodes, and then move this particle to one of them */

      for(i = 0, n = 0; i<2; i++)
	for(j=0; j<2; j++)
	  for(k=0; k<2; k++)
	    {
	      sun = get_empty_node();
	      current->suns[n++] = sun;
	      sun->len = 0.5 * current->len;
	      sun->center[0] = current->center[0] + 0.25 * (2*i-1) * current->len;
	      sun->center[1] = current->center[1] + 0.25 * (2*j-1) * current->len;
	      sun->center[2] = current->center[2] + 0.25 * (2*k-1) * current->len;
	    }

      /* determine in which subnode the old particle sits */
      p_subnode = get_subnode_index(current, current->p);
      
      /* now move the particle to this subnode */
      current->suns[p_subnode]->p = current->p;
      current->p = NULL;

      /* determine in which subnode the new particle sits */
      pnew_subnode = get_subnode_index(current, pnew);
      
      /* now try to insert the new particle there */
      insert_particle(current->suns[pnew_subnode], pnew);
    }
  else
    {
      /* check in which subnode the new particle would fall */
      pnew_subnode = get_subnode_index(current, pnew);
      
      /* if the corresponding subnode exists, we try to insert the particle there,
	 otherwise we know there are no subnodes in the node, so we can put the particle into the current node */
      if(current->suns[pnew_subnode])
	insert_particle(current->suns[pnew_subnode], pnew);
      else
	current->p = pnew;
    }
}


/* This function recursively calculates the multipole moments for the current node.
 * We only use monopoles here.
 */
void calc_multipole_moments(node *current)
{
  int n, j;

  if(current->suns[0])   /* do we have subnodes? */
    {
      /* yes, so let's first calculate their multipole moments */
      for(n = 0; n < 8; n++)
	calc_multipole_moments(current->suns[n]);

      /* initialize the node multipole moments to zero */
      current->mass  = 0;
      for(j = 0; j < 3; j++)
	current->cm[j] = 0;

      /* now calculate the moment of the current node from those of its suns */
      /*
       * ..... TO BE FILLED IN ....
       */
    }    
  else
    {
      if(current->p)  /* do we at least have a particle? */
	{
	  /* yes, so let's copy this particle to the multipole moments of the node */

	  current->mass = current->p->mass;
	  for(j = 0; j < 3; j++)
	    current->cm[j] = current->p->pos[j];
	}
      else
	{
	  /* nothing in here at all; let's initialize the multipole moments to zero */
	  current->mass  = 0;
	  for(j = 0; j < 3; j++)
	    current->cm[j] = 0;
	}
    }
}


double get_opening_angle(node *current, double pos[3])
{
  int j;
  double r2 = 0;

  for(j = 0; j < 3; j++)
    r2 += (current->cm[j] - pos[j]) * (current->cm[j] - pos[j]);

  return current->len / (sqrt(r2) + 1.0e-35);
}



void walk_tree(node *current, double pos[3], double acc[3])
{
  int n;
  double theta;
  
  if(current->mass)   /* only do something if there is mass in this branch of the tree (i.e. if it is not empty) */
    {
      theta = get_opening_angle(current, pos);
  
      /* if the node is seen under a small enough angle or contains a single particle,
       * we take its multipole expansion, and we're done for this branch 
       */
      if(theta < opening_threshold || current->p) 
	{
	  /*
	   * ..... TO BE FILLED IN ....
	   *
	   *     acc[0] += ....
	   *     acc[1] += ....
	   *     acc[2] += ....
	   */
	}
      else
	{
	  /* otherwise we have to open the node and look at all daughter nodes in turn */

	  if(current->suns[0])             /* make sure that we actually have subnodes */
	    for(n=0; n<8; n++)
	      walk_tree(current->suns[n], pos, acc);
	}
    }
}



int main(int argc, char **argv)
{
  node *root;
  int i, j, N;
  double t0, t1;

  N = MAX_POINTS;

  srand48(42);   /* set a random number seed */

  /* create a random particle set, uniformly distributed in a box */
  for(i=0; i < N; i++)
    {
      star[i].mass = 1.0 / N;

      for(j=0; j<3; j++)
	star[i].pos[j] = drand48();      
    }

  /* create an empty root node for the tree */
  root = get_empty_node();

  /* set the dimension and position of the root node */
  root->len = 1.0;
  for(j=0; j<3; j++)
    root->center[j] = 0.5;
 
  /* insert the particles into the tree */
  for(i=0; i < N; i++)
    insert_particle(root, &star[i]);

  /* calculate the multipole moments */
  calc_multipole_moments(root);


  /* set a timer */
  t0 = (double) clock();

  /* now calculate the accelerations with the tree */
  for(i = 0; i < N; i++)
    {
      for(j = 0; j < 3; j++)
	star[i].acc_tree[j] = 0;

      walk_tree(root, star[i].pos, star[i].acc_tree);
    }

  t1 = (double) clock();
  printf("\nforce calculation with tree took:        %8g sec\n", (t1 - t0) / CLOCKS_PER_SEC);


  t0 = (double) clock();

  /* now calculate the accelerations with direct summation, for comparison */
  for(i = 0; i < N; i++)
    {
	  /*
	   * ..... TO BE FILLED IN ....
	   *
	   *     star[i].acc_exact[0] = ....
	   *     star[i].acc_exact[0] = ....
	   *     star[i].acc_exact[0] = ....
	   */
    }
  
  t1 = (double) clock();
  printf("\ncalculation with direct summation took:  %8g sec\n", (t1 - t0) / CLOCKS_PER_SEC);

  /* now do the calculation of the mean relative error 
   */

  double err_sum = 0;

  /*
   * ..... TO BE FILLED IN ....
   *
   */

  err_sum /= N;

  printf("\nAverage relative error:  %8g\n", err_sum);

  exit(0);
}
