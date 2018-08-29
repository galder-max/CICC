void hardass(int *N, int *ass, int *mat)
{
  int n=N[0];
  for(int i=0;i<n;i++)
      for(int j=0;j<n;j++)
	  if(ass[j]==ass[i])
	    mat[i*n+j]=1;
}
