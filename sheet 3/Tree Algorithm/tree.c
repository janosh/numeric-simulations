#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <time.h>


// Let's define two types of structures, one for the particles, one for the tree nodes
typedef struct particle_data {
    double pos[3];
    double mass;

    double acc_tree[3];
    double acc_exact[3];
} particle;

typedef struct node_data {
    double center[3];
    double length;

    double com[3];
    double mass;
    double quadrupole[3][3];

    struct node_data *subnds[8];

    particle *p;
} node;


// A static global variable or a function has scope only in the file it's declared in. A static variable inside a function keeps its value between invocations.
static const double eps = 0.001;    // This is the gravitational softening length.
static const double G = 1.0;        // Gravitational cosntant
static int count_nodes = 0;         // Enumerates the active nodes during a given calculation. It is reset to zero before each change of the calculation parameters.


// This function returns a pointer to an empty tree node.
node *get_empty_node(node *tree, int num_nodes)
{
    if(count_nodes < num_nodes)
    {
        node *node_ptr = &tree[count_nodes++];
        memset(node_ptr, 0, sizeof(node));    // The C library function void *memset(void *str, int c, size_t n) copies the character c (an unsigned char) to the first n characters of the string pointed to by the argument str.
        return node_ptr;
    }
    else
    {
        printf("Sorry, we are out of tree nodes. Aborting calculation.\n");
        exit(1);
    }
}


// This function determines the index (0-7) of the subnode in which a particle falls within a given node.
int get_subnode_index(node *current, particle *p)
{
    int index = 0;

    if(p->pos[0] > current->center[0])  // p->pos gets the member called pos from the struct that p points to.
    index += 4;
    if(p->pos[1] > current->center[1])
    index += 2;
    if(p->pos[2] > current->center[2])
    index += 1;

    return index;
}


// This function's task is to insert a new particle into a given tree node.
void insert_particle(node *current, particle *pnew, node *tree, int num_nodes)
{
    node *subnd;
    int p_subnode, pnew_subnode;

    if(current->p)  // Does the node already contain a particle? If so, we need to create a new set of 8 subnodes, and then move this particle to one of them.

    {
        for(int i = 0, n = 0; i<2; i++)
        for(int j=0; j<2; j++)
        for(int k=0; k<2; k++)
        {
            subnd = get_empty_node(tree, num_nodes);
            current->subnds[n++] = subnd;   // This assigns the above newly created subnode to the n-th entry of the struct array subnds[8] of type node. subnds[8] belongs to the parent node current which already contained a particle.
            subnd->length = 0.5 * current->length;
            subnd->center[0] = current->center[0] + 0.25 * (2*i-1) * current->length;
            subnd->center[1] = current->center[1] + 0.25 * (2*j-1) * current->length;
            subnd->center[2] = current->center[2] + 0.25 * (2*k-1) * current->length;
        }

        // Here, we determine in which subnode the old particle was.
        p_subnode = get_subnode_index(current, current->p);

        // We then move the particle to this subnode.
        current->subnds[p_subnode]->p = current->p;
        current->p = NULL;  // NULL is a macro for the null pointer.

        // Next, we determine which subnode the new particle occupies.
        pnew_subnode = get_subnode_index(current, pnew);

        // Finally, we try to insert the new particle there.
        insert_particle(current->subnds[pnew_subnode], pnew, tree, num_nodes);
    }
    else
    {
        // The current node was empty. We now check into which subnode the new particle will fall.
        pnew_subnode = get_subnode_index(current, pnew);

        // If the corresponding subnode exists, we try to insert the particle there. If it does not, we know there are no subnodes in the current node, so we can put the particle into the current node directly.
        if(current->subnds[pnew_subnode])
        insert_particle(current->subnds[pnew_subnode], pnew, tree, num_nodes);
        else
        current->p = pnew;
    }
}


double get_dist_to_com(node *current, int k) {
    double dist_to_com;
    for (int i = 0; i < 3; i++) {
        dist_to_com += (current->com[i] - current->subnds[k]->p->pos[i]) * (current->com[i] - current->subnds[k]->p->pos[i]);
    }
    return dist_to_com;
}


// This function recursively calculates the multipole moments for the current node. We only use monopoles here.
void calc_multipole_moments(node *current)
{
    if(current->subnds[0])   // Do we have subnodes?
    {
        // Yes, so let's first calculate their multipole moments
        for(int n = 0; n < 8; n++)
        calc_multipole_moments(current->subnds[n]);

        // We initialize the node multipole moments to zero
        current->mass  = 0;
        for(int j = 0; j < 3; j++)
        current->com[j] = 0;

        // Now, we can calculate the moment of the current node from those of its subnodes.
        // Monopole:
        for (int k = 0; k < 8; k++)
        current->mass += current->subnds[k]->mass;

        // Next, we calculate the center of mass.
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < 8; k++) {
                current->com[i] += current->subnds[k]->mass * current->subnds[k]->com[i];
            }
            current->com[i] /= current->mass;
        }
        //
    }
    else    // There are no subnodes.
    {
        if(current->p)  // Do we at least have a particle?
        {
            // Yes, so we copy this particle to the multipole moments of the node.

            current->mass = current->p->mass;
            for(int j = 0; j < 3; j++)
            current->com[j] = current->p->pos[j];
        }
        else
        {
            // There is nothing in this node at all so we initialize the multipole moments to zero.
            current->mass  = 0;
            for(int j = 0; j < 3; j++)
            current->com[j] = 0;
        }
    }
}


double get_opening_angle(node *current, double pos[3])
{
    double r2 = 0;
    for(int j = 0; j < 3; j++)
    r2 += (current->com[j] - pos[j]) * (current->com[j] - pos[j]);

    return current->length / (sqrt(r2) + 1.0e-35);  // We add here the double-at-zero machine epsilon in the denominator to avoid 'nan' when r2 = 0.
}



void walk_tree(node *current, double pos[3], double acc[3], int *interaction_count, double opening_angle)
{
    double theta;

    if(current->mass)   // Only do something if there is mass in this branch of the tree (i.e. if it is not empty).
    {
        theta = get_opening_angle(current, pos);

        // If the node is seen under a small enough angle or contains a single particle, we take its multipole expansion, and we're done for this branch.
        if(theta < opening_angle || current->p)
        {
            double r2 = 0;
            for(int j = 0; j < 3; j++) {
                r2 += (pos[j] - current->com[j]) * (pos[j] - current->com[j]);
            }
            for (int i = 0; i < 3; i++) {
                acc[i] += -G * current->mass * (pos[i] - current->com[i]) / pow(r2 + eps * eps, 1.5);
            }
            // Here, we count the total number of particle-node interactions.
            (*interaction_count)++;
        }
        else    // Otherwise we have to open the node and look at all daughter nodes in turn.
        {
            if(current->subnds[0])             // We have to make sure that there actually are any subnodes. Technically, the case of no subnodes and no particle should not arise for any node.
            for(int n=0; n<8; n++)
            walk_tree(current->subnds[n], pos, acc, interaction_count, opening_angle);
        }
    }
}



int main(int argc, char **argv)
{
    // srand48(42);   // This sets a random number seed. More precisely, it initializes the uniformly distributed double-precision pseudo-random number generator rand48(). The seed is used to create a random particle set, uniformly distributed in a box.
    srand48(clock());   // To avoid getting the same pesudo-random number sequence every time, initialize the sequence with a varying source by passing srand the system time. This requires <time.h>.

    FILE *N_data = fopen("N.dat", "w");
    FILE *tree_data = fopen("tree.dat", "w");
    FILE *direct_data = fopen("direct.dat", "w");
    FILE *error_data = fopen("error.dat", "w");
    FILE *interaction_data = fopen("interaction.dat", "w");

    for (int a = 0; a < 4; a++) {   // Loop over different particle numbers N.
        for (int b = 0; b < 3; b++) {   // Loop over different opening angles.
            int N = 5000 * pow(2, a);   // This sets the number of particles to be used in the current calculation.
            int num_nodes = 5 * N;   // This sets the number of nodes to be used in the current calculation.
            double opening_angle = 0.2 * pow(2, b);      // We define the tree opening angle.


            // Let's create, for simplicity, some static arrays which will hold our data.
            node        tree[num_nodes];
            particle    star[N];
            printf("\nCalculation %d.%d\nN = \t\t%d\ntheta = \t%g°\n", a+1, b+1, N, opening_angle);
            if (b == 1) {
                fwrite(&N, sizeof(int), 1, N_data);
            }
            for(int i=0; i < N; i++)
            {
                star[i].mass = 1.0 / N; // Total mass M normalized to 1. Every particle gets the same mass of m = 1/N.

                for(int j=0; j<3; j++)
                star[i].pos[j] = drand48(); // drand48() uses a linear congruential algorithm and 48-bit integer arithmetic.
                // Note here that the dot operator and arrow operator in C are very similarly. (*star).pos is the same as star->pos. However, the dot operator can't be applied to pointers.
            }

            // This creates an empty root node for the tree.
            node *root = get_empty_node(tree, num_nodes);

            // We then set the dimension and position of the root node.
            root->length = 1.0; // We assign the root node a unit length.
            for(int j=0; j<3; j++)
            root->center[j] = 0.5; // And position the node's center right in the middle.

            // The next step is to insert the particles into the tree.
            for(int i=0; i < N; i++)
            insert_particle(root, &star[i], tree, num_nodes);

            // Once the particles exist, we are able to calculate the multipole moments.
            calc_multipole_moments(root);


            // This starts a timer for the tree method calculation.
            double t0 = (double) clock();

            // We now calculate the accelerations with the tree.
            int interaction_count = 0;
            for(int i = 0; i < N; i++)
            {
                for(int j = 0; j < 3; j++)
                star[i].acc_tree[j] = 0;

                walk_tree(root, star[i].pos, star[i].acc_tree, &interaction_count, opening_angle);
            }
            // This ends the timer started by t0.
            double t1 = (double) clock();
            double tree_exec_time = (t1 - t0) / CLOCKS_PER_SEC;
            printf("t_tree = \t%g sec\n", tree_exec_time);
            if (b == 1) {
                fwrite(&tree_exec_time, sizeof(double), 1, tree_data);
            }

            // This starts a timer for the direct summation calculation.
            t0 = (double) clock();

            // For comparison with our tree method, we calculate the accelerations via direct summation.
            for(int i=0; i<N; i++) {
                for(int l=0; l<3; l++) {
                    star[i].acc_exact[l]=0;
                }
            }

            for(int i = 0; i < N; i++) {
                for(int j = 0; j < N; j++) {
                    double r2 = 0;
                    for(int k = 0; k < 3; k++) {
                        r2 += (star[i].pos[k] - star[j].pos[k]) * (star[i].pos[k] - star[j].pos[k]);
                    }
                    for (int l = 0; l < 3; l++) {
                        star[i].acc_exact[l] += -G * star[j].mass * (star[i].pos[l] - star[j].pos[l]) / pow(r2 + eps * eps, 1.5);
                    }
                }
            }

            t1 = (double) clock();
            double direct_exec_time = (t1 - t0) / CLOCKS_PER_SEC;
            printf("t_dir = \t%g sec\n", direct_exec_time);
            if (b == 1) {
                fwrite(&direct_exec_time, sizeof(double), 1, direct_data);
            }
            // Finally, we wish to calculate the mean relative error.

            double err_sum = 0, acc_diff=0, acc_denom=0;

            // We calculate the total error for all particles here.
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < 3; j++) {
                    acc_diff += (star[i].acc_tree[j] - star[i].acc_exact[j]) * (star[i].acc_tree[j] - star[i].acc_exact[j]);
                    acc_denom += star[i].acc_exact[j] * star[i].acc_exact[j];
                }
                err_sum += sqrt(acc_diff)/sqrt(acc_denom);
            }
            // And then reduce it to a per-particle average.
            double avrg_error = err_sum/N, avrg_int_count = (double) interaction_count/N;
            printf("eta = \t\t%g\nint_count = \t%g\n\n", avrg_error, avrg_int_count);
            if (b == 1) {
                fwrite(&avrg_error, sizeof(double), 1, error_data);
                fwrite(&avrg_int_count, sizeof(double), 1, interaction_data);
            }

            // We perform a spot check comparison of forces obtained via tree method and direct summation (just one particle) to verify that the force approximation is roughly correct.
            int spot_check = (int) N * drand48();
            printf("Spot check for particle %d\n", spot_check);
            for (int i = 0; i < 3; i++) {
                printf("a_tree[%d] = \t%.6g\na_direct[%d] = \t%.6g\n", i+1, star[spot_check].acc_tree[i], i+1, star[spot_check].acc_exact[i]);
            }

            count_nodes = 0;    // This resets the node counter to zero to begin a fresh count in the next calculcation.
        }   // End of for loop over different opening angles.
    }   // End of for loop over different particle numbers.



    exit(0);
}
