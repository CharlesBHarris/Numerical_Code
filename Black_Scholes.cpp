//
//  main.cpp
//  
//
//  Created by charles harris on 9/21/14.
//  Feel free to steal it.
//



///////////////////////////////////////////////////////////////////////////
//                                                                       //
//           Black-Scholes Implicit Untransformed dx-dt                  //
//           American Put BC, non-dividend.                              //
//           See p.456 of Options, Futures, and Other Derivatives 8th    //
//           by John Hull for details.                                   //         
//                                                                       //
///////////////////////////////////////////////////////////////////////////



#include <iostream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <sstream>




using namespace std;
int main(int argc, const char * argv[])
{
    
    
    
    
    //declare variables
    int N=10, M=20;
    double S[N+1][M+1];
    double K=50.0,Smax=100.0,T=.4166666666667,r=0.1;
    double dt=T/N,ds=Smax/M;
    double a[M-1],b[M-1],c[M-1],d[M-1],p[M],q[M];
    double sigma=.4;
    
    
    
    //set boundary values
    for (int i=0;i<=M;i++)
    {
        if (K-i*ds>0)
            S[N][i]=K-i*ds;
        
    }
    for (int i=0;i<=N;i++)
    {
        S[i][0]=K;
        S[i][M]=0;
    }
    
    //calculate coefficients
    for (int i=0;i<=M-2;i++)
    {
        a[i]=0.5*r*(i+1)*dt-0.5*sigma*sigma*(i+1)*(i+1)*dt;
        b[i]=1+sigma*sigma*(i+1)*(i+1)*dt+r*dt;
        c[i]=-0.5*r*(i+1)*dt-0.5*sigma*sigma*(i+1)*(i+1)*dt;
    }
    
    //initialize d
    for (int i=0;i<=M-2;i++)
    {
        d[i]=S[N][i+1];
    }
    
    //calculate tridiagonal coefficients
    p[0]=0;
    q[0]=S[N-1][0];
    for (int i=0;i<=M-2;i++)
    {
        p[i+1]=-c[i]/(b[i]+a[i]*p[i]);
        q[i+1]=(d[i]-a[i]*q[i])/(b[i]+a[i]*p[i]);
    }
    
    //solve via Thomas algorithm
    for(int l=N-1;l>=0;l--)
    {
        for (int i=M-1;i>=1;i--)
        {
            S[l][i]=q[i]+p[i]*S[l][i+1];
        }
        //option loop
        for (int i=1;i<=M-1;i++)
        {
            if(S[l][i]<K-i*ds)
                S[l][i]=K-i*ds;
        }
    
        //reset d
        for (int i=0;i<=M-2;i++)
        {
            d[i]=S[l][i+1];
        }
        
        //reset q
        q[0]=S[l-1][0];
        for (int i=0;i<=M-2;i++)
        {
            q[i+1]=(d[i]-a[i]*q[i])/(b[i]+a[i]*p[i]);
        }
        
    }
    
    //output
    for(int i=0;i<=M;i++)
    {
        cout<< S[0][i] <<endl;
    }
    
       return 0;
}





