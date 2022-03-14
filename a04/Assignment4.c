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

#define MINPROB      1.0e-12
#define SOURCE0      0
#define SOURCE1      2
#define DELTA		 2
#define ENDTIME		 3600
#define PI           3.1415926

typedef struct tproxel *pproxel;

typedef struct tproxel {
    int     id;                  /* unique proxel id for searching    */
    int     s;                   /* discrete state of SPN             */
    int     tau1k;               /* first supplementary variable      */
    int     tau2k;               /* second supplementary variable     */
    double  val;                 /* proxel probability                */
    pproxel left, right;         /* pointers to child proxels in tree */
} proxel;

double  *y[3];             /* vectors for storing solution      */
double  tmax       ;         /* maximum simulation time           */
int     TAUMAX;
int     totcnt;                 /* counts total proxels processed    */
int     maxccp;                 /* counts max # concurrent proxels   */
int     ccpcnt;                 /* counts concurrent proxels         */
proxel *root[2];                /* trees for organising proxels      */
proxel *firstfree  = NULL;      /* linked list of free proxels       */
double  eerror      = 0;         /* accumulated error                 */
int     sw         = 0;         /* switch for old and new time steps */
int len;
double dt;
int cntsrc0;                    /* count for source0 proxels */
int cntsrc1;                    /* count for source1 proxels */

/********************************************************/
/*	distribution functions			                    */
/*	instantaneous rate functions			            */
/********************************************************/
 

/* returns weibull IRF */
double weibullhrf(double x, double alpha, double beta, double x0) {
    double y;

    y = beta/alpha*pow((x-x0)/alpha,beta-1);
    return(y);
}


/* returns deterministic IRF */
double dethrf(double x, double d) {
    double y;

    if (fabs(x - d) < dt/2)
        y = 1.0/dt;
    else
        y = 0.0;
    return(y);

}


/* returns uniform IRF */
double unihrf(double x, double a, double b) {
    double y;

    if ((x >= a) && (x < b))
        y = 1.0/(b-x);
    else
        y = 0.0;

    return(y);
}

/* returns exponential IRF */
double exphrf(double x, double l) {
    return(l);
}

double normalpdf(double x, double m, double s){
   double z = (x-m)/s;

   return (exp(-z * z/2)/(sqrt(2*PI)*s));
}


double logGamma(double x){
   double coef[] = {76.18009173, -86.50532033, 24.01409822, -1.231739516, 0.00120858003, -0.00000536382};
   double step = 2.50662827465, fpf = 5.5, t, tmp, ser;
   int i;

   t = x - 1;
   tmp = t + fpf;
   tmp = (t + 0.5) * log(tmp) -tmp;
   ser = 1;
   for(i = 1; i <= 6; i++){
      t = t+1;
      ser = ser + coef[i-1]/t;
   }
   return (tmp+log(step * ser));
}

double gammaSeries(double x, double a){
   int n, maxit = 100;
   double eps = 0.0000003;
   double sum = 1.0 / a, ap = a, gln = logGamma(a), del = sum;

   for (n = 1; n <= maxit; n++){
      ap++;
      del = del * x / ap;
      sum = sum + del;
      if (fabs(del) < fabs(sum) * eps) break;
   }
   return (sum * exp(-x + a * log(x) - gln));
}

double gammaCF(double x, double a){
   int n, maxit = 100;
   double eps = 0.0000003;
   double gln = logGamma(a), g = 0, gold = 0, a0 = 1, a1 = x, b0 = 0, b1 = 1, fac = 1;
   double an, ana, anf;

   for(n = 1; n <= maxit; n++){
      an = 1.0 * n;
      ana = an - a;
      a0 = (a1 + a0 * ana) * fac;
      b0 = (b1 + b0 * ana) * fac;
      anf = an * fac;
      a1 = x * a0 + anf * a1;
      b1 = x * b0 + anf * b1;
      if (a1 != 0){
         fac = 1.0 / a1;
         g = b1 * fac;
         if (fabs((g-gold)/g) < eps) 
            break;
         gold = g;
      }
   }
   return (exp(-x + a * log(x) - gln) * g);
}

double gammacdf(double x, double a){
   if (x <= 0) 
      return 0;
   else 
   if (x < a + 1) 
      return gammaSeries(x,a);
   else 
      return (1 - gammaCF(x, a));
}

double normalcdf(double x, double m, double s){
   double z = (x - m) / s;

   if (z >= 0)
      return 0.5 + 0.5 * gammacdf(z * z / 2, 0.5);
   else 
      return (0.5 - 0.5 * gammacdf(z * z / 2, 0.5));
}

/* returns normal IRF */
double normalhrf(double x, double m, double s){
   return(normalpdf(x, m, s)/(1 - normalcdf(x, m, s)));
}

double lognormalpdf(double x, double a, double b){
   double z = (log(x) - a) / b;

   return(exp(- z * z / 2) / (x * sqrt(2 * PI) * b));
}

double lognormalcdf(double x, double a, double b){
    double z = (log(x) - a) / b;

    if(x == 0)
        return 0;
    if (z >= 0) 
        return(0.5 + 0.5 * gammacdf(z * z / 2, 0.5));
    else 
        return(0.5 - 0.5 * gammacdf(z * z / 2, 0.5));
}

/* returns lognormal IRF using mu & sigma */
double lognormalhrf(double x, double a, double b){
   if ((x == 0.0) || (x > 70000))
      return(0.0);
   else
      return(lognormalpdf(x, a, b) / (1.0 - lognormalcdf(x, a, b)));
}

/********************************************************/
/*	output functions			                        */
/********************************************************/

/* print all proxels in tree */
void printtree(proxel *p) {
    if (p == NULL)
        return;
    printf("s %d t1 %d t2 %d val %lf \n",p->s,p->tau1k,p->tau2k,p->val);
    printtree(p->left);
    printtree(p->right);
}

/* print out complete solution */
void plotsolution(int kmax) {
	printf("\n\n");
    int k;

    for(k=0; k<=kmax; k++)
        printf("%7.5lf\t%7.5le\t%7.5le\n", k*dt, y[0][k], y[2][k]);

}

/* print out a proxel */
void printproxel(proxel *c) {
    printf("processing %2d %2d %7.5le \n", c->s, c->tau1k, c->val);
}

/********************************************************/
/*	proxel manipulation functions			            */
/********************************************************/

/* compute unique id from proxel state */
int state2id(int s, int t1k, int t2k) {
    return(TAUMAX*(TAUMAX*s+t1k)+t2k);
}

/* compute size of tree */
int size(proxel *p) {
    int sl, sr;

    if (p == NULL)
        return(0);
    sl = size(p->left);
    sr = size(p->right);
    return(sl+sr+1);
}

/* returns a proxel from the tree */
proxel *getproxel()
{
    proxel *temp;
    proxel *old;
    int LEFT = 0, RIGHT = 1;
    int dir, cont = 1;

    if (root[1-sw] == NULL)
        return(NULL);
    temp = root[1-sw];
    old  = temp;

    /* move down the tree to a leaf */
    while (cont == 1)
    {
        /* go right */
        if ((temp->right != NULL) && (temp->left == NULL))
        {
            old  = temp;
            temp = temp->right;
            dir  = RIGHT;
        }
        /* go left */
        else if ((temp->right == NULL) && (temp->left != NULL))
        {
            old  = temp;
            temp = temp->left;
            dir  = LEFT;
        }
        /* choose right/left at random */
        else if ((temp->right != NULL) && (temp->left != NULL))
        {
            if (rand() > RAND_MAX/2)
            {
                old  = temp;
                temp = temp->left;
                dir  = LEFT;
            }
            else
            {
                old  = temp;
                temp = temp->right;
                dir  = RIGHT;
            }
        }
        else
            cont = 0;
    }
    if (temp == root[1-sw])
        root[1-sw] = NULL;
    else
    {
        if (dir == RIGHT)
            old->right = NULL;
        else
            old->left  = NULL;
    }
    old = firstfree;
    firstfree = temp;
    temp->right = old;
    ccpcnt -= 1;
    return(temp);
}

/* get a fresh proxel and copy data into it */
proxel *insertproxel(int s, int tau1k, int tau2k, double val) {
    proxel *temp;

    /* create new proxel or grab one from free list */
    if (firstfree == NULL)
        temp = malloc(sizeof(proxel));
    else {
        temp = firstfree;
        firstfree = firstfree->right;
    }
    /* copy values */
    temp->id    = state2id(s, tau1k, tau2k);
    temp->s     = s;
    temp->tau1k = tau1k;
    temp->tau2k = tau2k;
    temp->val   = val;
    ccpcnt     += 1;
    if (maxccp < ccpcnt) {
        maxccp = ccpcnt;
        //printf("\n ccpcnt=%d",ccpcnt);
    }
    return(temp);
}

/* adds a new proxel to the tree */
void addproxel(int s, int tau1k, int tau2k, double val) {
    proxel *temp, *temp2;
    int cont = 1,id;

    /* Alarm! TAUMAX overstepped! */
    if (tau1k >= TAUMAX) {
        //  printf(">>> %3d %3d %3d %7.5le \n", s, tau1k, val, TAUMAX);
        tau1k = TAUMAX - 1;
    }


    /* New tree, add root */
    if (root[sw] == NULL) {
        root[sw] = insertproxel(s,tau1k, tau2k, val);
        root[sw]->left = NULL;
        root[sw]->right = NULL;
        return;
    }

    /* compute id of new proxel */
    id = state2id(s,tau1k, tau2k);

    /* Locate insertion point in tree */
    temp = root[sw];
    while (cont == 1) {
        if ((temp->left != NULL) && (id < temp->id))
            temp = temp->left;
        else
            if ((temp->right != NULL) && (id > temp->id))
                temp = temp->right;
            else
                cont = 0;
    }

    /* Insert left leaf into tree */
    if ((temp->left == NULL) && (id < temp->id)) {
        temp2        = insertproxel(s, tau1k,tau2k, val);
        temp->left   = temp2;
        temp2->left  = NULL;
        temp2->right = NULL;
        return;
    }

    /* Insert right leaf into tree */
    if ((temp->right == NULL) && (id > temp->id)) {
        temp2        = insertproxel(s, tau1k,tau2k, val);
        temp->right  = temp2;
        temp2->left  = NULL;
        temp2->right = NULL;
        return;
    }

    /* Proxels have the same id, just add their vals */
    if (id == temp->id) {
        temp->val += val;
        return;
    }
    printf("\n\n\n!!!!!! addproxel failed !!!!!\n\n\n");
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
/*  main processing loop                                */
/********************************************************/

int main(int argc, char **argv) {
    printf("===========================timestep = %d============================\n", DELTA);
    int     k, j, kmax;
    proxel *currproxel;
    double  val, z;
    int     s, tau1k; //, tau2k;

    /* initialise the simulation */
    root[0] = NULL;
    root[1] = NULL;
    eerror=0.0;
    totcnt  = 0;
    maxccp  = 0;
	double tmax = ENDTIME;
    dt = DELTA;
    kmax=tmax/dt+1;
    for (k = 0; k < 3; k++) {
        y[k] = malloc(sizeof(double)*(kmax+2));
        for (j = 0; j < kmax+2; j++) 
        	y[k][j] = 0.0;
    }
    TAUMAX = tmax/dt+1;
 
    /* set initial proxel */
    addproxel(1, 0, 0, 1.0);
    sw = 1 - sw;
    totcnt++;
    currproxel = getproxel();
    val        = currproxel->val;
    tau1k      = currproxel->tau1k;
    /*tau2k      = currproxel->tau2k;*/
    z = dt * sourceZero(tau1k * dt);
    addproxel(SOURCE0, 0, 0, get_log(val) + get_log(z));
    z = dt * sourceOne(tau1k * dt);
    addproxel(SOURCE1, 0, 0, get_log(val) + get_log(z));

    
    /* first loop: iteration over all time steps*/
    for (k = 2; k <= kmax; k++) {
        
        //if(k==79 || k==78) printtree(root[sw]);
        
         //printf("\nSTEP %d\n",k);
        /* current model time is k*dt */
        
        /* print progress information */
        if (k%100==0)  {
            //printf("\nSTEP %d\n",k);
            //printf("Size of tree %d\n",size(root[sw]));
        }
        
        sw = 1 - sw;

        /* second loop: iterating over all proxels of a time step */
        while (root[1-sw] != NULL)
        {
            totcnt++;
            currproxel = getproxel();
            while ((currproxel->val < MINPROB) && (root[1-sw] != NULL)) {
                val=currproxel->val;
                eerror += val;
                currproxel = getproxel();
            }
            val        = currproxel->val;
            tau1k      = currproxel->tau1k;
            /*tau2k      = currproxel->tau2k;*/
            s          = currproxel->s;
            y[s][k-1] += val;
            
            /* create child proxels */
            switch (s) {
            case SOURCE0:
                cntsrc0++;
                z = dt * sourceZero(tau1k * dt);
                addproxel(SOURCE0, tau1k + dt, 0, get_log(val) + get_log(z));
                z = dt * sourceOne(tau1k * dt);
                addproxel(SOURCE1, 0, 0, get_log(val) + get_log(z));
                break;
            case SOURCE1:
                cntsrc1++;
                z = dt * sourceZero(tau1k * dt);
                addproxel(SOURCE0, 0, 0, get_log(val) + get_log(z));
                z = dt * sourceOne(tau1k * dt);
                addproxel(SOURCE1, tau1k + dt, 0, get_log(val) + get_log(z));
                break;
            }
        }
    }
    printf("error = %7.5le\n", eerror);
    printf("ccpx = %d\n", maxccp);
    printf("count = %d\n", totcnt);
    printf("count of source 0 = %d\n", cntsrc0);
    printf("count of source 1 = %d\n", cntsrc1);
    printf("average defected parts by source 0  = %f\n", cntsrc0*0.1);
    printf("average defected parts by source 1  = %f\n", cntsrc1*0.05);
    printf("====================================================================\n");
    //plotsolution(kmax);
    return(0);
}