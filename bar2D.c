#include <stdio.h>
#include <stdlib.h>
struct element
{
    float length;
    float area;
    float stiffness;
    float bodyforcevalue;
    float tractionforcevalue;
    int globalnode_i;
    int globalnode_j;
    float kmatrix[2][2];
    float bodyforce[2][1];
    float tractionforce[2][1];
    float stress[1][1];
    float strain[1][1];
    float lmatrix[1][2];
    float qmatrix[2][1];
};
int main()
{
    int noele;
    printf("enter the number of elements present\n");
    scanf("%d",&noele);
    struct element elementdata[noele];
    int i=0;
    int j=0;
    int k,x,y,z;
    for(i==0;i<noele;i++)
    {
        printf("enter the datas of element %d\n",i+1);
        printf("enter the length of element %d\n",i+1);
        scanf("%f",&elementdata[i].length);
        printf("enter the area of element %d\n",i+1);
        scanf("%f",&elementdata[i].area);
        printf("enter the stiffness of element %d\n",i+1);
        scanf("%f",&elementdata[i].stiffness);
        printf("enter the global node number for element %d\n",i+1);
        printf("enter the global node i value for element%d\n",i+1);
        scanf("%d",&elementdata[i].globalnode_i);
        printf("enter the global node j value for element%d\n",i+1);
        scanf("%d",&elementdata[i].globalnode_j);
        printf("enter the body force acting of element %d, enter 0 if no force is acting\n",i+1);
        scanf("%f",&elementdata[i].bodyforcevalue);
        printf("enter the traction force acting of element %d, enter 0 if no force is acting\n",i+1);
        scanf("%f",&elementdata[i].tractionforcevalue);
    }
    int no_of_pt_force=0;
    printf("enter the number of external point forces present\n");
    scanf("%d",&no_of_pt_force);
    float pointforcevector[no_of_pt_force][2];
    for(i=0;i<no_of_pt_force;i++)
    {
        printf("enter the global node number at which point force is acting \n");
        scanf("%f",&pointforcevector[i][0]);
        printf("enter the force value\n");
        scanf("%f",&pointforcevector[i][1]);
    }
    i=0;
    for(i=0;i<noele;i++)
    {
        elementdata[i].kmatrix[0][0]=elementdata[i].stiffness*elementdata[i].area/elementdata[i].length;
        elementdata[i].kmatrix[0][1]=-elementdata[i].stiffness*elementdata[i].area/elementdata[i].length;;
        elementdata[i].kmatrix[1][0]=-elementdata[i].stiffness*elementdata[i].area/elementdata[i].length;;
        elementdata[i].kmatrix[1][1]=elementdata[i].stiffness*elementdata[i].area/elementdata[i].length;;
        elementdata[i].bodyforce[0][0]=0.5*elementdata[i].area*elementdata[i].length*elementdata[i].bodyforcevalue;
        elementdata[i].bodyforce[1][0]=0.5*elementdata[i].area*elementdata[i].length*elementdata[i].bodyforcevalue;
        elementdata[i].tractionforce[0][0]=0.5*elementdata[i].length*elementdata[i].tractionforcevalue;
        elementdata[i].tractionforce[1][0]=0.5*elementdata[i].length*elementdata[i].tractionforcevalue;
    }
    i=0;
    int noofnode=0;
    for(i=0;i<noele;i++)
    {
        if(elementdata[i].globalnode_i>=noofnode)
        {
            noofnode=elementdata[i].globalnode_i;
        }
        if(elementdata[i].globalnode_j>=noofnode)
        {
            noofnode=elementdata[i].globalnode_j;
        }
    }

    float globalkmatrix[noofnode][noofnode];
    float globalforcevector[noofnode][1];
    float globaldisplacementvector[noofnode][1];
    i=0;
    for(i=0;i<noofnode;i++)
    {
        globalforcevector[i][0]=0;
    }
    i=0;
    for(i=0;i<noofnode;i++)
    {
        for(j=0;j<noofnode;j++)
        {
            globalkmatrix[i][j]=0;
        }
        j=0;
    }
    i=0;
    int aux=0;
    for(i==0;i<noele;i++)
    {
        globalkmatrix[elementdata[i].globalnode_i-1][elementdata[i].globalnode_i-1]=globalkmatrix[elementdata[i].globalnode_i-1][elementdata[i].globalnode_i-1]+elementdata[i].kmatrix[0][0];
        globalkmatrix[elementdata[i].globalnode_i-1][elementdata[i].globalnode_j-1]=globalkmatrix[elementdata[i].globalnode_i-1][elementdata[i].globalnode_j-1]+elementdata[i].kmatrix[0][1];
        globalkmatrix[elementdata[i].globalnode_j-1][elementdata[i].globalnode_i-1]=globalkmatrix[elementdata[i].globalnode_j-1][elementdata[i].globalnode_i-1]+elementdata[i].kmatrix[1][0];
        globalkmatrix[elementdata[i].globalnode_j-1][elementdata[i].globalnode_j-1]=globalkmatrix[elementdata[i].globalnode_j-1][elementdata[i].globalnode_j-1]+elementdata[i].kmatrix[1][1];
        //globalforcevector[elementdata[i].globalnode_i-1][0]=globalforcevector[elementdata[i].globalnode_i-1][0]+elementdata[i].tractionforce[0][0]+elementdata[i].bodyforce[0][0];
        //globalforcevector[elementdata[i].globalnode_j-1][0]=globalforcevector[elementdata[i].globalnode_j-1][0]+elementdata[i].tractionforce[1][0]+elementdata[i].bodyforce[1][0];
    }
    for(i=0;i<no_of_pt_force;i++)
    {
        aux=pointforcevector[i][0]-1;
        globalforcevector[aux][0]=globalforcevector[aux][0]+pointforcevector[i][1];
    }
    i=0;
    j=0;
    int bctype=0;
    float maximum=-100;
        i=0;
        for(i==0;i<noofnode;i++)
        {
            for(j==0;j<noofnode;j++)
            {
                if(maximum<=globalkmatrix[i][j])
                    maximum=globalkmatrix[i][j];
            }
            j=0;
        }
    float cmagnitude;
    cmagnitude=maximum*10000;
        int noofbc=0;
        printf("enter the number of boundary condition present\n");
        scanf("%d",&noofbc);
        float boundaryvalues[noofbc][2];
        for(i=0;i<noofbc;i++)
        {
            printf("enter the global node number at which boundary condition is specified\n");
            scanf("%f",&boundaryvalues[i][0]);
            printf("enter the displacement boundary condition\n");
            scanf("%f",&boundaryvalues[i][1]);
        }
        i=0;
        int dummy=0;
        for(i=0;i<noofbc;i++)
        {
            dummy=boundaryvalues[i][0]-1;
            globalkmatrix[dummy][dummy]=globalkmatrix[dummy][dummy]+cmagnitude;
            globalforcevector[dummy][0]=globalforcevector[dummy][0]+cmagnitude*boundaryvalues[i][1];
        }

    i=0;
    j=0;
        int noofmbc=0;
        printf("enter the number of multipoint constraint present,enter 0 if no multipoint constraint is there\n");
        scanf("%d",&noofmbc);
        float multipoint[noofmbc][5];
        for(i=0;i<noofmbc;i++)
        {
            printf("enter the first global node number at which multipoint boundary condition is specified\n");
            scanf("%f",&multipoint[i][0]);
            printf("enter the second global node number at which multipoint boundary condition is specified\n");
            scanf("%f",&multipoint[i][1]);
            printf("enter the beta_1 value\n");
            scanf("%f",&multipoint[i][2]);
            printf("enter the beta_2 value\n");
            scanf("%f",&multipoint[i][3]);
            printf("enter the beta_0 value\n");
            scanf("%f",&multipoint[i][4]);
        }
        int dummy1=0;
        int dummy2=0;
        for(i=0;i<noofmbc;i++)
        {
            dummy1=multipoint[i][0]-1;
            dummy2=multipoint[i][2]-1;
            globalkmatrix[dummy1][dummy1]=globalkmatrix[dummy1][dummy1]+cmagnitude*multipoint[i][2]*multipoint[i][2];
            globalkmatrix[dummy2][dummy2]=globalkmatrix[dummy2][dummy2]+cmagnitude*multipoint[i][3]*multipoint[i][3];
            globalkmatrix[dummy1][dummy2]=globalkmatrix[dummy1][dummy2]+cmagnitude*multipoint[i][2]*multipoint[i][3];
            globalkmatrix[dummy2][dummy1]=globalkmatrix[dummy2][dummy1]+cmagnitude*multipoint[i][2]*multipoint[i][3];
            globalforcevector[dummy1][0]=globalforcevector[dummy1][0]+cmagnitude*multipoint[i][4]*multipoint[i][2];
            globalforcevector[dummy2][0]=globalforcevector[dummy2][0]+cmagnitude*multipoint[i][4]*multipoint[i][3];
        }

    float aug[noofnode][noofnode+1];
    for(x=0;x<(noofnode);x++)
        for(y=0;y<(noofnode);y++)
            aug[x][y]=globalkmatrix[x][y];
    for(x=0;x<(noofnode);x++)
        aug[x][noofnode]=globalforcevector[x][0];
    printf("augmented\n");
        for(x=0;x<noofnode;x++)
        {   for(y=0;y<noofnode+1;y++)
                printf("%f\t",aug[x][y]);
            printf("\n");
            printf("\n");

        }
    float c,sum;
    float globaldisp[noofnode];
    for(j=0; j<noofnode; j++)
    {
        for(i=0; i<noofnode; i++)
        {
            if(i>j)
            {
                c=aug[i][j]/aug[j][j];
                for(k=0; k<noofnode+1; k++)
                {
                    aug[i][k]=aug[i][k]-c*aug[j][k];
                }
            }
        }
    }


    globaldisp[noofnode-1]=(aug[noofnode-1][noofnode])/(aug[noofnode-1][noofnode-1]);
    for(i=noofnode-2; i>=0; i--)
    {
        sum=0;
        for(j=i+1; j<=noofnode-1; j++)
        {
            sum=sum+aug[i][j]*globaldisp[j];
        }
        globaldisp[i]=(aug[i][noofnode]-sum)/aug[i][i];
    }
    printf("\nThe solution is: \n");
    for(i=0; i<noofnode; i++)
    {
        printf("\nx%d=%f\t",i,globaldisp[i]);
    }

    float disp[noofnode][1];
    for(x=0;x<(noofnode);x++)
        disp[x][0]=globaldisp[x];
    for (i=0;i<noele;i++)
    {
         elementdata[i].qmatrix[0][0]=globaldisp[elementdata[i].globalnode_i-1];
         elementdata[i].qmatrix[1][0]=globaldisp[elementdata[i].globalnode_j-1];
    }
    for (i=0;i<noele;i++)
    {
        elementdata[i].lmatrix[0][0]=-1/(elementdata[i].length);
        elementdata[i].lmatrix[0][1]=1/(elementdata[i].length);

        sum=0;
        for (x=0;x<1;x++){
            for (y=0;y<1;y++){
                for (z=0;z<2;z++){
                    sum=sum + elementdata[i].lmatrix[x][z]*elementdata[i].qmatrix[z][y];
                }
            elementdata[i].strain[x][y]=sum;
            elementdata[i].stress[x][y]=elementdata[i].strain[x][y]*elementdata[i].stiffness;
            sum=0;
        }
      }
    }
    for (i=0;i<noele;i++)
    {
        printf("the stress in element %d is:\n",i+1);
        printf("%f\n",elementdata[i].stress[0][0]);
        printf("the strain in element %d is:\n",i+1);
        printf("%f\n",elementdata[i].strain[0][0]);
    }
    float kunmodifiedmatrix[noofnode][noofnode];
    for(x=0;x<(noofnode);x++)
        for(y=0;y<(noofnode);y++)
            kunmodifiedmatrix[x][y]=globalkmatrix[x][y];
    for(i=0;i<noofbc;i++)
        {
            dummy=boundaryvalues[i][0]-1;
            kunmodifiedmatrix[dummy][dummy]=kunmodifiedmatrix[dummy][dummy]-cmagnitude;
        }
    for(i=0;i<noofmbc;i++)
        {
            dummy1=multipoint[i][0]-1;
            dummy2=multipoint[i][2]-1;
            kunmodifiedmatrix[dummy1][dummy1]=kunmodifiedmatrix[dummy1][dummy1]-cmagnitude*multipoint[i][2]*multipoint[i][2];
            kunmodifiedmatrix[dummy2][dummy2]=kunmodifiedmatrix[dummy2][dummy2]-cmagnitude*multipoint[i][3]*multipoint[i][3];
            kunmodifiedmatrix[dummy1][dummy2]=kunmodifiedmatrix[dummy1][dummy2]-cmagnitude*multipoint[i][2]*multipoint[i][3];
            kunmodifiedmatrix[dummy2][dummy1]=kunmodifiedmatrix[dummy2][dummy1]-cmagnitude*multipoint[i][2]*multipoint[i][3];
        }
    sum=0;
    float reactionforcevector[noofnode][1];
        for (x=0;x<noofnode;x++){
            for (y=0;y<1;y++){
                for (z=0;z<noofnode;z++){
                    sum=sum + kunmodifiedmatrix[x][z]*disp[z][y];
                }
            reactionforcevector[x][y]=sum;
            sum=0;
            }

        }
        for(x=0;x<(noofnode);x++)
            reactionforcevector[x][0]=reactionforcevector[x][0]-globalforcevector[x][0];

    printf("reaction force is\n");
        for(x=0;x<noofnode;x++)
        {   for(y=0;y<1;y++)
                printf("%f\t",reactionforcevector[x][y]);
            printf("\n");
            printf("\n");

        }
    return 0;
}
