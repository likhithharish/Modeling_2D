 Attachments
  •  Scanned by Gmail
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define pi 3.14159265

int n,netnode,store[100],l,store1[100],select;
double Kfinal[100][100],fmodified[100][1],displacement[100][1],Kinverse[100][100],di[100];


struct frame {
int node1,node2,dof1[3],dof2[3];
double area,moi,x,y,x1,y1,angle,length,ym,height;
double stiffness[6][6],Lt[6][6],L[6][6];
double moment,stress,strain,shear,dof3[6][1];
}key[100];

struct nodes{
double k1,k2,k3,k4;
 int nde,a[3],b[3];
}q[100],q1[100];


double find(double x,double y,double x1,double y1)
{
   double change,val;
    val= 180.0/pi;

  if((x!=x1)&&(y!=y1))
    {
     if(((x<x1)&&(y<y1))||((x>x1)&&(y>y1)))
       change=atan(abs(y-y1)/abs(x-x1))*val;
      else
       change=180-(atan(abs(y-y1)/abs(x-x1))*val);
    }
    else if(x==x1)
        change=90;
    else
        change=0;

   return change;
}


void element_def()
{
 int i,count=0;
 printf("Define the elements:\n");

  for(i=1;i<=n;i++)
  {
      printf("\nElement:%d",i);
      printf("\ncoordinates(x,y)and(x1,y1):\n");
      scanf("%lf%lf%lf%lf",&key[i].x,&key[i].y,&key[i].x1,&key[i].y1);
      ++count;
      q[count].k1=key[i].x;
      q[count].k2=key[i].y;
      ++count;
      q[count].k1=key[i].x1;
      q[count].k2=key[i].y1;
      printf("Length:\n");
      scanf("%lf",&key[i].length);
      printf("Height:\n");
      scanf("%lf",&key[i].height);
      printf("Area of cross-section and Moment of inertia:\n");
      scanf("%lf%lf",&key[i].area,&key[i].moi);
      printf("Elasticity modulus:\n");
      scanf("%lf",&key[i].ym);
      key[i].angle=find(key[i].x,key[i].y,key[i].x1,key[i].y1);
  }
}

void classify()
{
  int flag=0,flag1=1,i,j;

  for(i=1;i<=(2*n);i++)
   {
    if (i==1)
    {
    ++flag;
    q1[flag].k3=q[i].k1;
    q1[flag].k4=q[i].k2;
    q1[flag].nde=flag;
    }
   else
   {
    for(j=1;j<=flag;j++)
     {
      if((q[i].k1==q1[j].k3)&&(q[i].k2==q1[j].k4))
        flag1=0;
     }
    if(flag1==1)
    {
      ++flag;
      q1[flag].k3=q[i].k1;
      q1[flag].k4=q[i].k2;
      q1[flag].nde=flag;
    }
  }
  flag1=1;
   }

 netnode=flag;

 int k=0;
 for(i=1;i<=flag;i++)
 {
  for(j=1;j<=3;j++)
    q1[i].a[j]= ++k;
 }

 int h;
 for(i=1;i<=n;i++)
 {
     for(j=1;j<=flag;j++)
     {
         if((key[i].x==q1[j].k3)&&(key[i].y==q1[j].k4))
         {key[i].node1=q1[j].nde;
         for(h=1;h<=3;h++)
            key[i].dof1[h]=q1[j].a[h];
         }
     }
     for(j=1;j<=flag;j++)
     {
        if((key[i].x1==q1[j].k3)&&(key[i].y1==q1[j].k4)){
           key[i].node2=q1[j].nde;
         for(h=1;h<=3;h++)
            key[i].dof2[h]=q1[j].a[h];
          }
     }
 }


}

void multiply(int no,double lte[6][6],double ki[6][6],double lt[6][6])
{
  double multiply[6][6],sum=0;
  int i,j,k;

    for (i=0;i<6;i++){
     for (j=0;j<6;j++){
      for (k=0;k<6;k++){
        sum=sum+(lte[i][k]*ki[k][j]);
      }
       multiply[i][j]=sum;
       sum=0;
      }
    }

  for (i=0;i<6;i++){
     for (j=0;j<6;j++){
      for (k=0;k<6;k++){
        sum=sum+(multiply[i][k]*lt[k][j]);
      }
       key[no].stiffness[i][j]=sum;
       sum=0;
     }
    }
}

double Kstiff[100][100];

void stiff (int i)
{
    double k1,k2,k3,k4,k5,a=12,b=6,c=4,d=2;
    int m,n,i1,j1;

    k1=(key[i].ym*key[i].area)/key[i].length;
    k2=a*((key[i].ym*key[i].moi)/pow(key[i].length,3));
    k3=b*((key[i].ym*key[i].moi)/pow(key[i].length,2));
    k4=c*((key[i].ym*key[i].moi)/key[i].length);
    k5=d*((key[i].ym*key[i].moi)/key[i].length);

    double theta1,theta2,theta;
    theta=(key[i].angle*pi)/180;
    theta1=cos(theta);
    theta2=sin(theta);
    double Ke[6][6]={{k1,0,0,-k1,0,0},{0,k2,k3,0,-k2,k3},{0,k3,k4,0,-k3,k5},{-k1,0,0,k1,0,0},{0,-k2,-k3,0,k2,-k3},{0,k3,k5,0,-k3,k4}};
    double L[6][6]={{theta1,theta2,0,0,0,0},{-theta2,theta1,0,0,0,0},{0,0,1,0,0,0},{0,0,0,theta1,theta2,0},{0,0,0,-theta2,theta1,0},{0,0,0,0,0,1}};
    double Lt1[6][6]={{theta1,-theta2,0,0,0,0},{theta2,theta1,0,0,0,0},{0,0,1,0,0,0},{0,0,0,theta1,-theta2,0},{0,0,0,theta2,theta1,0},{0,0,0,0,0,1}};



    for(i1=0;i1<6;i1++)
    {
        for(j1=0;j1<6;j1++)
        {
            key[i].Lt[i1][j1]=Lt1[i1][j1];
            key[i].L[i1][j1]=L[i1][j1];

    }
    }

    multiply(i,Lt1,Ke,L);

    Kstiff[(3*key[i].node1)-3][(3*key[i].node1)-3]=Kstiff[(3*key[i].node1)-3][(3*key[i].node1)-3]+key[i].stiffness[0][0];
    Kstiff[(3*key[i].node1)-3][(3*key[i].node1)-2]=Kstiff[(3*key[i].node1)-3][(3*key[i].node1)-2]+key[i].stiffness[0][1];
    Kstiff[(3*key[i].node1)-3][(3*key[i].node1)-1]=Kstiff[(3*key[i].node1)-3][(3*key[i].node1)-1]+key[i].stiffness[0][2];
    Kstiff[(3*key[i].node1)-2][(3*key[i].node1)-3]=Kstiff[(3*key[i].node1)-2][(3*key[i].node1)-3]+key[i].stiffness[1][0];
    Kstiff[(3*key[i].node1)-2][(3*key[i].node1)-2]=Kstiff[(3*key[i].node1)-2][(3*key[i].node1)-2]+key[i].stiffness[1][1];
    Kstiff[(3*key[i].node1)-2][(3*key[i].node1)-1]=Kstiff[(3*key[i].node1)-2][(3*key[i].node1)-1]+key[i].stiffness[1][2];
    Kstiff[(3*key[i].node1)-1][(3*key[i].node1)-3]=Kstiff[(3*key[i].node1)-1][(3*key[i].node1)-3]+key[i].stiffness[2][0];
    Kstiff[(3*key[i].node1)-1][(3*key[i].node1)-2]=Kstiff[(3*key[i].node1)-1][(3*key[i].node1)-2]+key[i].stiffness[2][1];
    Kstiff[(3*key[i].node1)-1][(3*key[i].node1)-1]=Kstiff[(3*key[i].node1)-1][(3*key[i].node1)-1]+key[i].stiffness[2][2];

    Kstiff[(3*key[i].node1)-3][(3*key[i].node2)-3]=Kstiff[(3*key[i].node1)-3][(3*key[i].node2)-3]+key[i].stiffness[0][3];
    Kstiff[(3*key[i].node1)-3][(3*key[i].node2)-2]=Kstiff[(3*key[i].node1)-3][(3*key[i].node2)-2]+key[i].stiffness[0][4];
    Kstiff[(3*key[i].node1)-3][(3*key[i].node2)-1]=Kstiff[(3*key[i].node1)-3][(3*key[i].node2)-1]+key[i].stiffness[0][5];
    Kstiff[(3*key[i].node1)-2][(3*key[i].node2)-3]=Kstiff[(3*key[i].node1)-2][(3*key[i].node2)-3]+key[i].stiffness[1][3];
    Kstiff[(3*key[i].node1)-2][(3*key[i].node2)-2]=Kstiff[(3*key[i].node1)-2][(3*key[i].node2)-2]+key[i].stiffness[1][4];
    Kstiff[(3*key[i].node1)-2][(3*key[i].node2)-1]=Kstiff[(3*key[i].node1)-2][(3*key[i].node2)-1]+key[i].stiffness[1][5];
    Kstiff[(3*key[i].node1)-1][(3*key[i].node2)-3]=Kstiff[(3*key[i].node1)-1][(3*key[i].node2)-3]+key[i].stiffness[2][3];
    Kstiff[(3*key[i].node1)-1][(3*key[i].node2)-2]=Kstiff[(3*key[i].node1)-1][(3*key[i].node2)-2]+key[i].stiffness[2][4];
    Kstiff[(3*key[i].node1)-1][(3*key[i].node2)-1]=Kstiff[(3*key[i].node1)-1][(3*key[i].node2)-1]+key[i].stiffness[2][5];

    Kstiff[(3*key[i].node2)-3][(3*key[i].node1)-3]=Kstiff[(3*key[i].node2)-3][(3*key[i].node1)-3]+key[i].stiffness[3][0];
    Kstiff[(3*key[i].node2)-3][(3*key[i].node1)-2]=Kstiff[(3*key[i].node2)-3][(3*key[i].node1)-2]+key[i].stiffness[3][1];
    Kstiff[(3*key[i].node2)-3][(3*key[i].node1)-1]=Kstiff[(3*key[i].node2)-3][(3*key[i].node1)-1]+key[i].stiffness[3][2];
    Kstiff[(3*key[i].node2)-2][(3*key[i].node1)-3]=Kstiff[(3*key[i].node2)-2][(3*key[i].node1)-3]+key[i].stiffness[4][0];
    Kstiff[(3*key[i].node2)-2][(3*key[i].node1)-2]=Kstiff[(3*key[i].node2)-2][(3*key[i].node1)-2]+key[i].stiffness[4][1];
    Kstiff[(3*key[i].node2)-2][(3*key[i].node1)-1]=Kstiff[(3*key[i].node2)-2][(3*key[i].node1)-1]+key[i].stiffness[4][2];
    Kstiff[(3*key[i].node2)-1][(3*key[i].node1)-3]=Kstiff[(3*key[i].node2)-1][(3*key[i].node1)-3]+key[i].stiffness[5][0];
    Kstiff[(3*key[i].node2)-1][(3*key[i].node1)-2]=Kstiff[(3*key[i].node2)-1][(3*key[i].node1)-2]+key[i].stiffness[5][1];
    Kstiff[(3*key[i].node2)-1][(3*key[i].node1)-1]=Kstiff[(3*key[i].node2)-1][(3*key[i].node1)-1]+key[i].stiffness[5][2];

    Kstiff[(3*key[i].node2)-3][(3*key[i].node2)-3]=Kstiff[(3*key[i].node2)-3][(3*key[i].node2)-3]+key[i].stiffness[3][3];
    Kstiff[(3*key[i].node2)-3][(3*key[i].node2)-2]=Kstiff[(3*key[i].node2)-3][(3*key[i].node2)-2]+key[i].stiffness[3][4];
    Kstiff[(3*key[i].node2)-3][(3*key[i].node2)-1]=Kstiff[(3*key[i].node2)-3][(3*key[i].node2)-1]+key[i].stiffness[3][5];
    Kstiff[(3*key[i].node2)-2][(3*key[i].node2)-3]=Kstiff[(3*key[i].node2)-2][(3*key[i].node2)-3]+key[i].stiffness[4][3];
    Kstiff[(3*key[i].node2)-2][(3*key[i].node2)-2]=Kstiff[(3*key[i].node2)-2][(3*key[i].node2)-2]+key[i].stiffness[4][4];
    Kstiff[(3*key[i].node2)-2][(3*key[i].node2)-1]=Kstiff[(3*key[i].node2)-2][(3*key[i].node2)-1]+key[i].stiffness[4][5];
    Kstiff[(3*key[i].node2)-1][(3*key[i].node2)-3]=Kstiff[(3*key[i].node2)-1][(3*key[i].node2)-3]+key[i].stiffness[5][3];
    Kstiff[(3*key[i].node2)-1][(3*key[i].node2)-2]=Kstiff[(3*key[i].node2)-1][(3*key[i].node2)-2]+key[i].stiffness[5][4];
    Kstiff[(3*key[i].node2)-1][(3*key[i].node2)-1]=Kstiff[(3*key[i].node2)-1][(3*key[i].node2)-1]+key[i].stiffness[5][5];

    printf("\nElement (%d):\n",i);
   for(m=0;m<6;m++)
   {
       for(n=0;n<6;n++)
       {
           printf("%lf\t",key[i].stiffness[m][n]);
       }
       printf("\n");
   }
printf("\n");
}

double force[100][1];

void load()
{

double p,a,b,c;
int tot,i,j,k,arr[100],sum=0,i1;

printf("\nApplication of pressure:\n");
printf("No.of.elements on which pressure is acting:\n");
scanf("%d",&tot);
if(tot!=0)
{
  printf("\nEnter the element no.s:\n");
  for(i=1;i<=tot;i++)
    scanf("%d",&arr[i]);

for(i=1;i<=tot;i++)
{
    printf("\nPressure on element(%d): ",arr[i]);
    scanf("%lf",&p);
    double fe[6][1]={{0},{(p*key[arr[i]].length)/2},{(p*pow(key[arr[i]].length,2))/12},{0},{(p*key[arr[i]].length)/2},{-1*((p*pow(key[arr[i]].length,2))/12)}};
    double fei[6][1];

    for (i1=0;i1<6;i1++){
      for (k=0;k<6;k++){
        sum=sum+(key[arr[i]].Lt[i1][k]*fe[k][0]);
      }
       fei[i1][0]=sum;
       sum=0;
      }

    force[(3*key[arr[i]].node1)-3][0]=force[(3*key[arr[i]].node1)-3][0]+fei[0][0];
    force[(3*key[arr[i]].node1)-2][0]=force[(3*key[arr[i]].node1)-2][0]+fei[1][0];
    force[(3*key[arr[i]].node1)-1][0]=force[(3*key[arr[i]].node1)-1][0]+fei[2][0];
    force[(3*key[arr[i]].node2)-3][0]=force[(3*key[arr[i]].node2)-3][0]+fei[3][0];
    force[(3*key[arr[i]].node2)-2][0]=force[(3*key[arr[i]].node2)-2][0]+fei[4][0];
    force[(3*key[arr[i]].node2)-1][0]=force[(3*key[arr[i]].node2)-1][0]+fei[5][0];
}
}

int arr1[100];

  printf("\nApplication of point Load/moment:\n");
  printf("\nNo.of.nodes on which point load or moment is acting:\n");
  scanf("%d",&tot);
  if(tot!=0)
{
  printf("\nEnter the node no.s:\n");
  for(i=1;i<=tot;i++)
    scanf("%d",&arr1[i]);

  for(i=1;i<=tot;i++)
  {
      printf("Node(%d):\n",arr1[i]);
      printf("Fx:");
      scanf("%lf",&a);
      force[(3*arr1[i])-3][0]=force[(3*arr[i])-3][0]+a;
      printf("\nFy:");
      scanf("%lf",&b);
      force[(3*arr1[i])-2][0]=force[(3*arr1[i])-2][0]+b;
      printf("\nMz:");
      scanf("%lf",&c);
      force[(3*arr1[i])-1][0]=force[(3*arr1[i])-1][0]+c;
  }
}


}

int gdof[100][1];
void boundary()
{
int i;
double a;
printf("No.of.degrees of freedom to be constrained:\n");
scanf("%d",&select);
if(select!=0)
{
    printf("\nEnter the d.o.f no.s:\n");
    for(i=0;i<select;i++)
        scanf("%d",&store1[i]);

    printf("\nEnter the displacement values:\n");

    for(i=0;i<select;i++)
    {
        printf("D.O.F(%d):",store1[i]);
        scanf("%lf",&di[i]);
    }

}
}


double determinant(double [][100], int);
void cofactor(double [][100],int ,double [][100]);
void transpose(double [][100],double [][100], int);

void secondary(int i)
{
    double point,zhi,k1,k2,k3,k4;
    printf("\nEnter the x-coordinate of element(%d) where values of secondary variables to be determined:\n",i);
    scanf("%lf",&point);
    zhi=((2*(point-key[i].x))/(key[i].x1-key[i].x))-1;
    k1=6*zhi;
    k2=(-2+(6*zhi))*(key[i].length/2);
    k3=-6*zhi;
    k4=(2+(6*zhi))*(key[i].length/2);

    double dh[1][4]={{k1,k2,k3,k4}};
    double dh1[1][4]={{6,6*(key[i].length/2),-6,6*(key[i].length/2)}};
    double g=0;

        g=g+(dh[0][0]*key[i].dof3[1][0]);
        g=g+(dh[0][1]*key[i].dof3[2][0]);
        g=g+(dh[0][2]*key[i].dof3[4][0]);
        g=g+(dh[0][3]*key[i].dof3[5][0]);


  key[i].moment=(4*key[i].ym*key[i].moi*g)/pow(key[i].length,2);

  key[i].shear=((6*key[i].ym*key[i].moi)*((2*key[i].dof3[1][0])+(key[i].length*key[i].dof3[2][0])-(2*key[i].dof3[4][0])+(key[i].length*key[i].dof3[5][0])))/pow(key[i].length,3);           //calculation of shear


double yp;
  printf("\nEnter the value of y for stress and strain determination:\n");
  scanf("%lf",&yp);

  double c;
  c=(key[i].ym/key[i].length)*(-key[i].dof3[0][0]+key[i].dof3[3][0]);
  key[i].stress=((-1*(key[i].moment*yp))/key[i].moi)+c;
  key[i].strain=key[i].stress/key[i].ym;


  printf("\nSecondary variable values for element(%d):\n",i);
  printf("\nMoment:%lf",key[i].moment);
  printf("\nShear stress:%lf",key[i].shear);
  printf("\nNormal stress:%lf",key[i].stress);
  printf("\nstrain:%lf",key[i].strain);

}

int main()
{
  int i,j,k;
  printf("No of elements in the system:\n");
  scanf("%d",&n);
  element_def();
  classify();

  printf("Total nodes and DOF:\n");

  for (i=0;i<(3*netnode);i++)
   {
     for (j=0;j<(3*netnode);j++)
      Kstiff[i][j]=0;
   }

  for(i=1;i<=n;i++)
  {
    printf("Element(%d):\n",i);
    printf("Nodes:%d,%d\n",key[i].node1,key[i].node2);
    printf("\nD.O.F(1):%d,%d,%d",key[i].dof1[1],key[i].dof1[2],key[i].dof1[3]);
    printf("\nD.O.F(2):%d,%d,%d\n",key[i].dof2[1],key[i].dof2[2],key[i].dof2[3]);
    printf("Angle:%lf\n",key[i].angle);
    stiff(i);
  }

  printf("Global stiffness matrix:\n");
  for(i=0;i<(3*netnode);i++)
  {
      for(j=0;j<(3*netnode);j++)
       {
          printf("%lf\t",Kstiff[i][j]);
       }
      printf("\n");

  }

  for(i=0;i<(3*netnode);i++)
     force[i][0]=0;
  load();                                           //

 printf("\nGlobal force vector:\n");
 for(i=0;i<(3*netnode);i++)
    printf("%lf\n",force[i][0]);

printf("\nBoundary conditions:\n");                        //
l=0;
int flag=0;

for(i=0;i<(3*netnode);i++)
{
    for(j=0;j<select;j++)
    {
      if((i+1)!=store1[j])
        ++flag;
    }

    if(flag==select)
       store[l++]=i;
    flag=0;
}

printf("\n Fixed D.O.F\n");
for(i=0;i<select;i++)
    printf("%d\t%lf\n",store1[i],di[i]);

printf("\nD.O.F to be determined:\n");
for(i=0;i<l;i++)
    printf("%d\n",store[i]+1);



for(i=0;i<l;i++)
{
 for(j=0;j<l;j++)
  Kfinal[i][j]=Kstiff[store[i]][store[j]];
}


 printf("\nModified stiffness matrix:\n");
  for(i=0;i<l;i++)
  {
      for(j=0;j<l;j++)
       {
          printf("%lf\t",Kfinal[i][j]);
       }
      printf("\n");

  }

 for(j=0;j<select;j++)
 {
 for(i=0;i<l;i++)
   fmodified[i][0]=force[store[i]][0]-((Kstiff[i][store1[j]-1])*di[j]);
 }


printf("\nModified force vector:\n");
 for(i=0;i<l;i++)
    printf("%lf\n",fmodified[i][0]);


 double b[100][100];
  double d;

  d = determinant(Kfinal,l);
  if (d == 0)
   printf("\nInverse of Entered Matrix is not possible\n");
  else
   cofactor(Kfinal,l,b);

   double sum=0;
   for(i=0;i<select;i++)
   {
    for(j=0;j<(3*netnode);j++)
       {
           if((store1[i]-1)==j)
            displacement[j][0]=di[i];
       }

   }

   for (i=0;i<l;i++){
      for (k=0;k<l;k++){
        sum=sum+(Kinverse[i][k]*fmodified[k][0]);
      }
       displacement[store[i]][0]=sum;
       sum=0;
    }

    for(i=1;i<=n;i++)                                             //
    {
        key[i].dof3[0][0]=displacement[(3*key[i].node1)-3][0];
        key[i].dof3[1][0]=displacement[(3*key[i].node1)-2][0];
        key[i].dof3[2][0]=displacement[(3*key[i].node1)-1][0];
        key[i].dof3[3][0]=displacement[(3*key[i].node2)-3][0];
        key[i].dof3[4][0]=displacement[(3*key[i].node2)-2][0];
        key[i].dof3[5][0]=displacement[(3*key[i].node2)-1][0];
    }

printf("\nD.O.F            Displacement-values:\n");
    for(i=0;i<(3*netnode);i++)
       printf("%d                    %lf\n",i+1,displacement[i][0]);

printf("\nSecondary variable calculations\n");
sum=0;

double multiply2[100][1];
 for (i=0;i<(3*netnode);i++){
    for (k=0;k<(3*netnode);k++){
        sum=sum+(Kstiff[i][k]*displacement[k][0]);
      }
       multiply2[i][0]=sum;
       sum=0;
    }

double reaction[100][1];
  for(i=0;i<(3*netnode);i++)
    reaction[i][0]=multiply2[i][0]-force[i][0];

printf("\nReaction-(forces and moments)\n");
printf("\nD.O.F            Reaction force/moments:\n");
    for(i=0;i<(3*netnode);i++)
       printf("%d              %lf\n",i+1,reaction[i][0]);

for(i=1;i<=n;i++)
    secondary(i);


return 0;
}

/*For calculating Determinant of the Matrix */
double determinant(double a[100][100],int k)
{
  double s = 1, det = 0,b[100][100];
  int i, j, m, n, c;
  if (k == 1)
    {
     return (a[0][0]);
    }
  else
    {
     det = 0;
     for (c = 0; c < k; c++)
       {
        m = 0;
        n = 0;
        for (i = 0;i < k; i++)
          {
            for (j = 0 ;j < k; j++)
              {
                b[i][j] = 0;
                if (i != 0 && j != c)
                 {
                   b[m][n] = a[i][j];
                   if (n < (k - 2))
                    n++;
                   else
                    {
                     n = 0;
                     m++;
                     }
                   }
               }
             }
          det = det + s * (a[0][c] * determinant(b, k - 1));
          s = -1 * s;
          }
    }

    return (det);
}

//Calculating cofactor
void cofactor(double num[100][100],int f,double b[100][100])
{
 double fac[100][100];
 int p, q, m, n, i, j;
 for (q = 0;q < f; q++)
 {
   for (p = 0;p < f; p++)
    {
     m = 0;
     n = 0;
     for (i = 0;i < f; i++)
     {
       for (j = 0;j < f; j++)
        {
          if (i != q && j != p)
          {
            b[m][n] = num[i][j];
            if (n < (f - 2))
             n++;
            else
             {
               n = 0;
               m++;
               }
            }
        }
      }
      fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
    }
  }
  transpose(num, fac, f);
}

/*Finding transpose of matrix*/
void transpose(double num[100][100], double fac[100][100], int r)
{
  int i, j;
  double b[100][100], d;

  for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
         b[i][j] = fac[j][i];
        }
    }
  d = determinant(num, r);
  for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
        Kinverse[i][j] = b[i][j] / d;
        }
    }
   printf("\n\n\nThe inverse of modified stiffness matrix is : \n");

   for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
         printf("\t%lf", Kinverse[i][j]);
        }
    printf("\n");
     }

}
