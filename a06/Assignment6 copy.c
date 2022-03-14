/*******************************************************************
 *                                                                 *
 *           Special Purpose Proxel-Based Solver                   *
 *                                                                 *
 *           Advanced Discrete Modelling 2006                      *
 *                                                                 *
 *           written/modified by                                   *
 *           Graham Horton, Sanja Lazarova-Molnar, Fabian Wickborn *
 ******************************************************************/

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include <sys/queue.h>
#include <string.h>
#include <stdbool.h>

#define MINPROB 1.0e-12
#define SOURCE0 0
#define SOURCE1 2
#define DELTA 10
#define ENDTIME 15050
#define PI 3.1415926

typedef struct tproxel *pproxel;

typedef struct tproxel
{
    int id;              /* unique proxel id for searching    */
    int s;               /* discrete state of SPN             */
    int tau1k;           /* first supplementary variable      */
    int tau2k;           /* second supplementary variable     */
    double val;          /* proxel probability                */
    pproxel left, right; /* pointers to child proxels in tree */
    char *path[1600];
    pproxel previous;
    int level;
    int active; /*for checking top prob proxels*/

} proxel;

double *y[3]; /* vectors for storing solution      */
double tmax;  /* maximum simulation time           */
int TAUMAX;
int totcnt;               /* counts total proxels processed    */
int maxccp;               /* counts max # concurrent proxels   */
int ccpcnt;               /* counts concurrent proxels         */
proxel *root;             /* trees for organising proxels      */
proxel *firstfree = NULL; /* linked list of free proxels       */
double eerror = 0;        /* accumulated error                 */
int sw = 0;               /* switch for old and new time steps */
int len;
double dt;
int curr_level = 0;

const int topN_count = 10; /* number of top N proxels to pick in each timestamp */
int index_topN = 0;         /* index for traversing proxels in each timestamp */
proxel *topN[20];          /* array for holding proxels in each timestamp */

/********************************************************/
/*	distribution functions			                    */
/*	instantaneous rate functions			            */
/********************************************************/

// A C program to demonstrate linked list based implementation of queue

// A linked list (LL) node to store a queue entry
struct QNode
{
    proxel *p;
    struct QNode *next;
};

// The queue, front stores the front node of LL and rear stores the
// last node of LL
struct Queue
{
    struct QNode *front, *rear;
};

// A utility function to create a new linked list node.
struct QNode *newNode(proxel *x)
{
    struct QNode *temp = (struct QNode *)malloc(sizeof(struct QNode));
    temp->p = x;
    temp->next = NULL;
    return temp;
}

// A utility function to create an empty queue
struct Queue *createQueue()
{
    struct Queue *q = (struct Queue *)malloc(sizeof(struct Queue));
    q->front = q->rear = NULL;
    return q;
}

// The function to add a key k to q
void enQueue(struct Queue *q, proxel *x)
{
    // Create a new LL node
    struct QNode *temp = newNode(x);

    // If queue is empty, then new node is front and rear both
    if (q->rear == NULL)
    {
        q->front = q->rear = temp;
        return;
    }

    // Add the new node at the end of queue and change rear
    q->rear->next = temp;
    q->rear = temp;
}

// Function to remove a key from given queue q
void deQueue(struct Queue *q)
{
    // If queue is empty, return NULL.
    if (q->front == NULL)
        return;

    // Store previous front and move front one node ahead
    struct QNode *temp = q->front;

    q->front = q->front->next;

    // If front becomes NULL, then change rear also as NULL
    if (q->front == NULL)
        q->rear = NULL;

    free(temp);
}

////////////////////////////////////

/* returns weibull IRF */
double weibullhrf(double x, double alpha, double beta, double x0)
{
    double y;

    y = beta / alpha * pow((x - x0) / alpha, beta - 1);
    return (y);
}

/* returns deterministic IRF */
double dethrf(double x, double d)
{
    double y;

    if (fabs(x - d) < dt / 2)
        y = 1.0 / dt;
    else
        y = 0.0;
    return (y);
}

/* returns uniform IRF */
double unihrf(double x, double a, double b)
{
    double y;

    if ((x >= a) && (x < b))
        y = 1.0 / (b - x);
    else
        y = 0.0;

    return (y);
}

/* returns exponential IRF */
double exphrf(double x, double l)
{
    return (l);
}

double normalpdf(double x, double m, double s)
{
    double z = (x - m) / s;

    return (exp(-z * z / 2) / (sqrt(2 * PI) * s));
}

double logGamma(double x)
{
    double coef[] = {76.18009173, -86.50532033, 24.01409822, -1.231739516, 0.00120858003, -0.00000536382};
    double step = 2.50662827465, fpf = 5.5, t, tmp, ser;
    int i;

    t = x - 1;
    tmp = t + fpf;
    tmp = (t + 0.5) * log(tmp) - tmp;
    ser = 1;
    for (i = 1; i <= 6; i++)
    {
        t = t + 1;
        ser = ser + coef[i - 1] / t;
    }
    return (tmp + log(step * ser));
}

double gammaSeries(double x, double a)
{
    int n, maxit = 100;
    double eps = 0.0000003;
    double sum = 1.0 / a, ap = a, gln = logGamma(a), del = sum;

    for (n = 1; n <= maxit; n++)
    {
        ap++;
        del = del * x / ap;
        sum = sum + del;
        if (fabs(del) < fabs(sum) * eps)
            break;
    }
    return (sum * exp(-x + a * log(x) - gln));
}

double gammaCF(double x, double a)
{
    int n, maxit = 100;
    double eps = 0.0000003;
    double gln = logGamma(a), g = 0, gold = 0, a0 = 1, a1 = x, b0 = 0, b1 = 1, fac = 1;
    double an, ana, anf;

    for (n = 1; n <= maxit; n++)
    {
        an = 1.0 * n;
        ana = an - a;
        a0 = (a1 + a0 * ana) * fac;
        b0 = (b1 + b0 * ana) * fac;
        anf = an * fac;
        a1 = x * a0 + anf * a1;
        b1 = x * b0 + anf * b1;
        if (a1 != 0)
        {
            fac = 1.0 / a1;
            g = b1 * fac;
            if (fabs((g - gold) / g) < eps)
                break;
            gold = g;
        }
    }
    return (exp(-x + a * log(x) - gln) * g);
}

double gammacdf(double x, double a)
{
    if (x <= 0)
        return 0;
    else if (x < a + 1)
        return gammaSeries(x, a);
    else
        return (1 - gammaCF(x, a));
}

double normalcdf(double x, double m, double s)
{
    double z = (x - m) / s;

    if (z >= 0)
        return 0.5 + 0.5 * gammacdf(z * z / 2, 0.5);
    else
        return (0.5 - 0.5 * gammacdf(z * z / 2, 0.5));
}

/* returns normal IRF */
double normalhrf(double x, double m, double s)
{
    return (normalpdf(x, m, s) / (1 - normalcdf(x, m, s)));
}

double lognormalpdf(double x, double a, double b)
{
    double z = (log(x) - a) / b;

    return (exp(-z * z / 2) / (x * sqrt(2 * PI) * b));
}

double lognormalcdf(double x, double a, double b)
{
    double z = (log(x) - a) / b;

    if (x == 0)
        return 0;
    if (z >= 0)
        return (0.5 + 0.5 * gammacdf(z * z / 2, 0.5));
    else
        return (0.5 - 0.5 * gammacdf(z * z / 2, 0.5));
}

/* returns lognormal IRF using mu & sigma */
double lognormalhrf(double x, double a, double b)
{
    if ((x == 0.0) || (x > 70000))
        return (0.0);
    else
        return (lognormalpdf(x, a, b) / (1.0 - lognormalcdf(x, a, b)));
}

/********************************************************/
/*	output functions			                        */
/********************************************************/

/* print all proxels in tree */
void printtree(proxel *p)
{
    if (p == NULL)
        return;
    printf("s %d t1 %d t2 %d val %lf \n", p->s, p->tau1k, p->tau2k, p->val);
    printtree(p->left);
    printtree(p->right);
}

/* print out complete solution */
void plotsolution(int kmax)
{
    printf("\n\n");
    int k;

    for (k = 0; k < kmax; k++)
        printf("%7.5lf\t%7.5le\t%7.5le\n", k * dt, y[0][k], y[2][k]);
}

/* print out a proxel */
void printproxel(proxel *c)
{
    printf("processing %2d %2d %7.5le \n", c->s, c->tau1k, c->val);
}

/********************************************************/
/*	proxel manipulation functions			            */
/********************************************************/

/* compute unique id from proxel state */
int state2id(int s, int t1k, int t2k)
{
    return (TAUMAX * (TAUMAX * s + t1k) + t2k);
}

/* compute size of tree */
int size(proxel *p)
{
    int sl, sr;

    if (p == NULL)
        return (0);
    sl = size(p->left);
    sr = size(p->right);
    return (sl + sr + 1);
}

/* returns a proxel from the tree */

/* get a fresh proxel and copy data into it */
proxel *insertproxel(int s, int tau1k, int tau2k, double val, int ct, pproxel prev)
{
    proxel *temp;

    /* create new proxel or grab one from free list */
    if (firstfree == NULL)
        temp = malloc(sizeof(proxel));
    else
    {
        temp = firstfree;
        firstfree = firstfree->right;
    }
    /* copy values */
    temp->previous = prev;
    temp->id = state2id(s, tau1k, tau2k);
    temp->s = s;
    temp->tau1k = tau1k;
    temp->tau2k = tau2k;
    temp->val = val;
    temp->level = ct;
    temp->active = 1;
    ccpcnt += 1;
    int i;
    if (temp->previous != NULL)
    {
        for (i = 0; i < ct; i++)
            temp->path[i] = temp->previous->path[i];

        temp->path[i] = s == 0 ? "SOURCE 0" : "SOURCE 1";
    }
    else
        temp->path[0] = s == 0 ? "SOURCE 0" : "SOURCE 1";

    if (maxccp < ccpcnt)
    {
        maxccp = ccpcnt;
        // printf("\n ccpcnt=%d",ccpcnt);
    }
    topN[index_topN++] = temp;
    return (temp);
}

/////NEW FUNCTIONS

/* adds a new proxel to the tree */
void addproxel(int s, int tau1k, int tau2k, double val, int curr_tstep, proxel *curr_proxel)
{

    if (root == NULL && curr_tstep == 0)
    {
        root = insertproxel(s, tau1k, tau2k, val, curr_tstep, NULL);
        root->left = NULL;
        root->right = NULL;
        root->previous = NULL;
        return;
    }

    proxel *temp2;
    int id = state2id(s, tau1k, tau2k);
    if ((curr_proxel->left == NULL))
    {
        temp2 = insertproxel(s, tau1k, tau2k, val, curr_tstep, curr_proxel);
        curr_proxel->left = temp2;
        temp2->left = NULL;
        temp2->right = NULL;
        temp2->previous = curr_proxel;
        return;
    }

    /* Insert right leaf into tree */
    if ((curr_proxel->right == NULL))
    {
        temp2 = insertproxel(s, tau1k, tau2k, val, curr_tstep, curr_proxel);
        curr_proxel->right = temp2;
        temp2->left = NULL;
        temp2->right = NULL;
        temp2->previous = curr_proxel;
        return;
    }
}

proxel *get(proxel *root)
{
    // If the tree is empty, assign new node address to root
    proxel *temp = malloc(sizeof(proxel));

    if (root->right == NULL && root->left == NULL)
        return root;

    // Else, do level order traversal until we find an empty
    // place, i.e. either left child or right child of some
    // node is pointing to NULL.

    struct Queue *q = createQueue();
    enQueue(q, root);

    // push(root);

    while (q->front != NULL)
    {
        proxel *temp = q->front->p;
        // pop();
        deQueue(q);

        if (temp->left != NULL)
            // push(temp->left);
            enQueue(q, temp->left);
        else if (temp->active == 1)
        {
            // temp->left = CreateNode(data);
            return temp;
        }

        if (temp->right != NULL)
            // push(temp->right);
            enQueue(q, temp->right);
        else if (temp->active == 1)
        {
            // temp->right = CreateNode(data);
            return temp;
        }
    }
}

double get_log(double probability)
{
    if (log(probability) == -INFINITY)
        return log(probability + MINPROB);
    else if (log(probability) == INFINITY || isnan(log(probability)))
        return 1.0;
    else
        return log(probability);
}

/********************************************************/
/*	model specific distribtuions	                    */
/********************************************************/

/* INSTANTANEOUS RATE FUNCTION 1 */
double sourceZero(double age)
{
    // return unihrf(age, 0.25, .5);
    // return normalhrf(age, 0.3, 0.1);
    return normalhrf(age, 150, 25);
}

/* INSTANTANEOUS RATE FUNCTION 2 */
double sourceOne(double age)
{
    // return exphrf(age, 2);
    // return dethrf(age, 0.5);
    return normalhrf(age, 120, 20);
}

/********************************************************/
/*  main processing loop                                */
/********************************************************/


int main(int argc, char **argv)
{

    double trace_prob[1600];
    FILE *file;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    int index = 0;

    FILE *fp;
    char ch;
    char ch_prev;

    fp = fopen(argv[0], "r");

    if (fp == NULL)
    {
        printf("File not found. \n");
    }
    else
    {
        ch = fgetc(fp);
        ch_prev = ch;
        printf("File is opening..... \n\n");
        while (ch != EOF)
        {

            if (ch == '\r')
            {
                if (ch_prev == 'K')
                    trace_prob[index] = 0.95;
                else if (ch_prev == 'D')
                    trace_prob[index] = 0.05;
                index++;
            }

            ch_prev = ch;
            ch = fgetc(fp);
        }
    }

    fclose(fp);

    int k, j, kmax;

    proxel *currproxel;
    double val, z;
    int s, tau1k, tau2k;

    /* initialise the simulation */
    root = NULL;
    // root[1] = NULL;
    eerror = 0.0;
    totcnt = 0;
    maxccp = 0;
    double tmax = ENDTIME;
    dt = DELTA;
    kmax = tmax / dt + 1;
    for (k = 0; k < 3; k++)
    {
        y[k] = malloc(sizeof(double) * (kmax + 2));
        for (j = 0; j < kmax + 2; j++)
            y[k][j] = 0.0;
    }
    TAUMAX = tmax / dt + 1;

    /* set initial proxel, Starting from SOURCE1 because it is mentioned as faster prodcution machine */
    addproxel(1, 0, 0, 1, 0, 0);
    currproxel = get(root);
    val = currproxel->val;
    tau1k = currproxel->tau1k;
    // tau2k      = currproxel->tau2k;
    s = currproxel->s;
    y[s][k - 1] += val;
    z = dt * sourceZero(tau1k * dt);
    addproxel(SOURCE0, 0, 0, val + get_log(z) + get_log(trace_prob[k]), 1, currproxel);
    z = dt * sourceOne(tau1k * dt);
    addproxel(SOURCE1, 0, 0, val + get_log(z) + get_log(trace_prob[k]), 1, currproxel);
                

    // int first = 0;
    /* first loop: iteration over all time steps*/
    for (k = 2; k <= kmax; k++)
    {

        // proxel * topN_proxels[5];
        // if(k==79 || k==78) printtree(root[sw]);

        // printf("\nSTEP %d\n",k);
        /* current model time is k*dt */

        /* print progress information */
        if (k % 100 == 0)
        {
            printf("\nSTEP %d\n", k);
            printf("Size of tree %d\n", size(root));
        }
        if (k == 2)
            currproxel = get(root);
        proxel *prev_currproxel = currproxel;
        index_topN = 0;
        // first = 0;
        // sw = 1 - sw;
        /* second loop: iterating over all proxels of a time step */
        while (currproxel->level == prev_currproxel->level)
        {
            if (currproxel->active == 0)
            {
                currproxel = get(root);
                continue;
            }

            prev_currproxel = currproxel;
            totcnt++;
            // first = 1;

            /*while ((currproxel->val < log(MINPROB)) && (root != NULL)) {
                val=currproxel->val;
                eerror += val;
                currproxel = get(root);
            }*/
            val = currproxel->val;
            tau1k = currproxel->tau1k;
            // tau2k      = currproxel->tau2k;
            s = currproxel->s;
            y[s][k - 1] += val;

            switch (s)
            {
            case SOURCE0:
                z = dt * sourceZero(tau1k * dt);
                addproxel(SOURCE0, tau1k + dt, 0, val + get_log(z) + get_log(trace_prob[k]), k, currproxel);
                z = dt * sourceOne(tau1k * dt);
                addproxel(SOURCE1, 0, 0, val + get_log(z) + get_log(trace_prob[k]), k, currproxel);
                break;
            case SOURCE1:
                z = dt * sourceZero(tau1k * dt);
                addproxel(SOURCE0, 0, 0, val + get_log(z) + get_log(trace_prob[k]), k, currproxel);
                z = dt * sourceOne(tau1k * dt);
                addproxel(SOURCE1, tau1k + dt, 0, val + get_log(z) + get_log(trace_prob[k]), k, currproxel);
                break;
            }

            currproxel = get(root);
            // printf("Timestep %d\n", k);
            // printf("Root level = %d\n", root->level);
        }
        if (pow(2, k) > topN_count && index_topN < (topN_count * 2))
        {
            for (int i = 0; i < index_topN - 1; i++)
            {
                for (int j = i + 1; j < index_topN; j++)
                {
                    if (topN[i]->val > topN[j]->val)
                    {
                        proxel *t = topN[i];
                        topN[i] = topN[j];
                        topN[j] = t;
                    }
                }
            }
            for (int i = 0; i < (pow(2, k) - topN_count); i++)
                topN[i]->active = 0;
        }
        else if (pow(2, k) > topN_count && index_topN == (topN_count * 2))
        {
            for (int i = 0; i < index_topN - 1; i++)
            {
                for (int j = i + 1; j < index_topN; j++)
                {
                    if (topN[i]->val > topN[j]->val)
                    {
                        proxel *t = topN[i];
                        topN[i] = topN[j];
                        topN[j] = t;
                    }
                }
            }
            for (int i = 0; i < topN_count; i++)
                topN[i]->active = 0;
        }

        // currproxel=currproxel->previous;
        // Front -=1 ;
        // printf("Timestep %d\n", k);
    }
    printf("Probability of output sequence: %f\n", topN[index_topN-1]->val);
    //printf("Lowest %f\n", topN[0]->val);

    for(int i=1; i<=10; i++)
    {
        printf("Path %d: ",i);
        for(int j=1; j<ENDTIME/DELTA; j++)
        {
            printf("%s -> ", topN[index_topN-i]->path[j]);
        }
        printf("END OF PATH \n");
    }
    // printf("Root level = %d\n", root->level);
    // printf("error = %7.5le\n", eerror);
    // printf("ccpx = %d\n", maxccp);
    // printf("count = %d\n", totcnt);
    //plotsolution(kmax);
    return (0);
}