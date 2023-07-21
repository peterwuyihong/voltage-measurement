#include<bits/stdc++.h>
#define ll long long
#define inf 1<<30
using namespace std;
const int MAXN=210;
const double eps=1e-12;
double a[MAXN][MAXN],x[MAXN];
bool frex[MAXN],vis[MAXN];
int Gauss(int equ,int var){
	for(int i=0;i<=var;i++) x[i]=0,frex[i]=1;
	int col=0,row;
	for(row=0;row<equ&&col<var;row++,col++){
		int mx_r=row;
		for(int i=row+1;i<equ;i++)
			if(abs(a[i][col])>abs(a[mx_r][col]))
				mx_r=i;
		if(mx_r!=row)
			for(int j=col;j<=var;j++)
				swap(a[row][j],a[mx_r][j]);
		if(fabs(a[row][col])<eps){row--;continue;}
		for(int i=row+1;i<equ;i++)
			if(fabs(a[i][col])>eps){
				double tmp=a[i][col]/a[row][col];
				for(int j=col;j<=var;j++)
					a[i][j]-=a[row][j]*tmp;
				a[i][col]=0;
			}
	}
	for(int i=row;i<equ;i++)
		if(fabs(a[i][col])>eps)
			return -1;
	if(var>row) return var-row;
	for(int i=var-1;i>=0;i--){
		double tmp=a[i][var];
		for(int j=i+1;j<var;j++)
			tmp-=x[j]*a[i][j];
		x[i]=tmp/a[i][i];
	}return 0;
}
int main(){
	int equ,var,n,m,k,Q,u,v;
	double w;
	scanf("%d%d%d%d",&n,&m,&k,&Q);
	equ=var=n+1;
	for(int i=0;i<k;i++){
		scanf("%d%lf",&u,&w);
		a[u][u]=1;a[u][var]=w;
		vis[u]=1;
	}
	a[0][0]=1;a[0][var]=0;vis[0]=1;
	for(int i=0;i<m;i++){
		scanf("%d%d%lf",&u,&v,&w);
		if(!vis[u]) a[u][u]+=1.0/w,a[u][v]-=1.0/w;
		if(!vis[v]) a[v][v]+=1.0/w,a[v][u]-=1.0/w;
	}
	for(int i=0;i<equ;i++)
		if(!vis[i]) a[i][var]=0;
	int res=Gauss(equ,var);
	while(Q--){
		scanf("%d%d",&u,&v);
		printf("%.12f\n",x[u]-x[v]);
	}
}