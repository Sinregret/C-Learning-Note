#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// 回收数组
void recycle(double ** arr);

// 复制数组副本
double ** copyArr(double ** arr,int row,int col);

// 将字符串转为二维数组
double ** transform(char * str,int * row_n,int * column_n);

// 打印二维数组
int printArr(double ** array,int row_n,int column_n);

// 行列式计算-代数余子式递归
double calc(double ** array,int n);

// calc - 简写
double calc2(double ** array,int n);

// 行列式计算 - 下三角
double calc3(double ** array,int n);

// 余子项（二维数组）
double ** cofactor(double ** arrIni,int n,int i,int j);

// 一维转二维
double ** upDim(double * arr,int row,int col);

// 矩阵转置
double ** transposition (double ** arr,int row,int col);

// 矩阵相乘
double ** matrixMulti(double ** a,int m,int s1,double ** b,int s2,int n);

// 伴随矩阵
double ** adjoint(double ** arr,int n);

// 伴随矩阵求逆矩阵
double ** inverseMatrix (double ** arr,int n);

// 初等变换逆矩阵
double ** inverseMatrix2 (double ** arr,int n);

// 矩阵 - 行最简
double ** simple(double ** array,int row,int col);

// 求矩阵的秩
int rank(double ** array,int row,int col);

// 入口函数
int main(int argc, char **argv)
{

    char * e = "1,2,4|\
				2,3,3|\
				2,3,5";
    int * row_n = (int *)malloc(4);
    int * column_n = (int *)malloc(4);
    double ** array = transform(e,row_n,column_n);
    
    /*
    char * e2 = "1,2,3|\
				6,5,6|\
				7,8,4";
    int * row_n2 = (int *)malloc(4);
    int * column_n2 = (int *)malloc(4);
    double ** array2 = transform(e2,row_n2,column_n2);
    */
    
    // 转置
    //double ** arr2 = transposition(array2,* row_n2,*column_n2 );
    //printArr(arr2,* column_n2,* row_n2);
    // 下三角行列式计算
    //double result = calc3(array,*row_n);
    //printf("### %f ###",result);
    // 矩阵相乘
    //matrixMulti(array,* row_n,* column_n,array2,* row_n2,* column_n2);
    // 伴随矩阵求逆矩阵
	double ** result = inverseMatrix2( array, * row_n);
	printArr(result,* row_n,* row_n);
	// 矩阵 - 行最简
	//double ** result = simple(array, * row_n,* column_n);
	//printArr(result,* row_n,* column_n);
	// 求矩阵的秩
	//printf("rank = %d\n",rank(array,* row_n,* column_n)) ;
    
    free(row_n);
    free(column_n);
    recycle(array);
    return 0;
}   

// 将字符串转为二维数组
double ** transform(char * str,int * row_n,int * column_n)
{
	if(* str == '\0') 
	{
		* row_n = 0;
		* column_n = 0;
		return NULL;
	}
	int c = 0;
	int l = 0;
	int column_num = 0;
	int row_num = 0;
	char * str1 = str;
	char * str2 = str;
	// 遍历，计算行列数
    do {
		if(* str1 == ',') c++;
		if(* str1 == '|') l++;
	} while(* ++str1 != '\0');
	row_num = l + 1;
	column_num = c/row_num + 1;
	* row_n = row_num;
    * column_n = column_num;
    
	// 构建二维指针（数组）
    double ** arr_xy = upDim(malloc(row_num*column_num*sizeof(double)),row_num,column_num);
    
	// 格式转换，顺序存储
	double * arr_x = arr_xy[0];
	* arr_x = atof(str2);
	while(* str2){
		if(* ++str2 != ',' && *str2 != '|') continue;
		* ++arr_x = atof(++str2);
	}
    
    return arr_xy;
}

// 打印二维数组
int printArr(double ** array,int row_n,int column_n)
{
	if(array == NULL) 
	{
		printf("\t[矩阵为空]\n");
		return 0;
	}
	
	printf("*-------------------------------------*\n");
	for(int i=0;i<row_n;i++){
		for(int j=0;j<column_n;j++)
		{
			if(j==0) printf("|");
			printf(" %f",array[i][j]);
			if(j<column_n-1) printf(",");
		}
		printf("|\n");
	}
    printf("*-------------------------------------*\n");
    return 0;
}

// 代数余字式 -递归
double calc(double ** array,int n)
{
	if(n==1) return array[0][0];
	if(n==2) return array[0][0]*array[1][1]-array[0][1]*array[1][0];
	
	double sum = 0;
	double ** temp = upDim(malloc(sizeof(double)*(n-1)*(n-1)),n-1,n-1);
	
	for(int i=0;i<n;i++)
	{
		// 保存余子项
		for(int x=1;x<n;x++)
		{
			int flag=0;
			for(int y=0;y<n;y++)
			{
				if(y==i)
				{
					 flag=1;
				}
				else
				{
					temp[x-1][y-flag]=array[x][y];
				}
			}
		}
		
		sum += pow(-1,i+2)*array[0][i]*calc(temp,n-1);
		// printArr(temp,n-1,n-1);
	}
	recycle(temp);
	return sum;
}

// calc 简化代码
double calc2(double ** array,int n)
{
	if(n==1) return array[0][0];
	if(n==2) return array[0][0]*array[1][1]-array[0][1]*array[1][0];
	
	double sum = 0;
	
	for(int i=0;i<n;i++)
	{
		// 固定计算第一行的代数余子式
		double ** temp = cofactor(array,n,0,i);
		sum += pow(-1,i+2)*array[0][i]*calc(temp,n-1);
		
		recycle(temp);
		// printArr(temp,n-1,n-1);
	}

	return sum;
}

// 行列式计算 - 下三角
double calc3(double ** arr,int n)
{
	// 数组副本
	double ** array = copyArr(arr,n,n);
	double result=1;
	int i =0;
	while(i < n)
	{
		/*
		 *  若 array[i][i]==0，
		 *	判断从i+1行开始，将第一个i列不为0的行，记为k行，k行加入第i行
		 * 	若 array[i+n][i]都为0，则结果为0
		**/
		if(array[i][i]==0)
		{
			for(int m=i+1;m<n;m++)
			{
				if(array[m][i]!=0)
				{
					for(int s=i;s<n;s++)
					{
						array[i][s]=array[i][s]+array[m][s];
					}
					break;
				}
			}
			if(array[i][i]==0) return 0;
		}
		
		// 化简为下三角
		for(int h=i+1;h<n;h++)
		{
			// 每行系数
			double r = -array[h][i]/array[i][i];
			for(int s=i;s<n;s++)
			{
				array[h][s]=array[h][s]+r*array[i][s];
			}
		}
		
		i++;
	}
	// 返回结果
	for(int s=0;s<n;s++)
	{
		result *= array[s][s];
	}
	
	// 回收副本
	recycle(array);
	return result;	
}

// 矩阵 - 行最简
double ** simple(double ** arr,int row,int col)
{
	// 数组副本
	double ** array = copyArr(arr,row,col);

	int i=0,j=0;
	while(i<row)
	{
		while(j<col)
		{
			/*
			 *  若 array[i][j]==0，
			 *	判断从i+1行开始，将第一个j列不为0的行，记为k行，k行加入第i行
			 * 	若 array[i+n][j]都为0，则 j++
			 *  循环，直到 array[i][j] 不为 0 或 j=col
			**/
			if(array[i][j]==0)
			{
				for(int m=i+1;m<row;m++)
				{
					if(array[m][j]!=0)
					{
						for(int s=j;s<col;s++)
						{
							array[i][s]=array[i][s]+array[m][s];
						}
						break;
					}
				}
			}
			if(array[i][j]!=0) {
				break;
			} else {
				j++;
			}
			
		}
		
		// a[i][j] 依然为零，则从i行开始，矩阵全为0
		if(j==col) break;
		
		if(array[i][j]!=1){
			for(int s=col-1;s>=j;s--)
			{
				array[i][s]=array[i][s]/array[i][j];
			}
		}
		
		// 将j列中，除 a[i][j] 全化为0
		for(int h=0;h<row;h++)
		{
			if(h==i || array[h][j]==0) continue;
			for(int s=col-1;s>=j;s--)
			{
				array[h][s]=array[h][s]-array[h][j]*array[i][s];
			}
		
		}
		i++;
		j++;
	}
	
	return array;
}

// 求矩阵的秩
int rank(double ** arr,int row,int col)
{
	// 数组副本
	double ** array = copyArr(arr,row,col);
	
	int r = 0;
	array = simple(array,row,col); 
	for(int i = 0;i<row;i++)
	{
		for(int j = 0;j<row;j++)
		{
			if(array[i][j]!=0)
			{
				r++;
				break;
			}
		}
	}
	
	// 回收副本
	recycle(array);
	
	return r;
}

// 一维数组（指针）转二维数组
double ** upDim(double * arr,int row,int col)
{
	double ** temp = malloc(sizeof(double *)*row);
	for(int i=0;i<row;i++)
    {
		temp[i]=arr+col*i;
	}
	return temp;
}


// 两矩阵相乘
double ** matrixMulti(double ** a,int m,int s1,double ** b,int s2,int n)
{
	if(m<1 || n<1 || s1<1 || s2<1) 
	{
		printf("Log[ERROR] - matrixMulti - 存在矩阵行列数为零\n");
		return NULL;
	}
	if(s1!=s2) 
	{
		printf("Log[ERROR] - matrixMulti - 左矩阵列数与右矩阵行数不相等 :%d,%d\n",s1,s2);
		return NULL;
	}
	// double * matrixCon = malloc(sizeof(double)*m*n);
	double ** matrixP = upDim(malloc(sizeof(double)*m*n),m,n);
	
	for(int i = 0;i<m;i++)
	{
		for(int j = 0;j<n;j++)
		{
			double sum =0;
			for(int k = 0;k<s1;k++){
				sum += a[i][k]*b[k][j];
			}
			matrixP[i][j] = sum;
		}
	}
	
	//printArr(matrixP,m,n);
	return matrixP;
}

// 矩阵转置
double ** transposition (double ** arr,int row,int col)
{
	double ** result = upDim(malloc(sizeof(double *)*row*col),row,col);
	for(int i = 0;i<row;i++)
	{
		for(int j = 0;j<col;j++)
		{
			result[i][j]=arr[j][i];
		}
	}
	return result;
}

// 伴随矩阵
double ** adjoint(double ** arr,int n)
{

	double ** reslut = upDim(malloc(sizeof(double)*n*n),n,n);
	
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			double ** tmp = cofactor(arr,n,i,j);
			reslut[j][i] = pow(-1,i+j)*calc2(tmp,n-1);
			recycle(tmp);
		}
	}
	return reslut;	
}

// arrIni[i][j]对应的余子项
double ** cofactor(double ** arrIni,int n,int i,int j)
{
	double ** result = upDim(malloc(sizeof(double)*(n-1)*(n-1)),n-1,n-1);
	int rowFlag=0;
	for(int x=0;x<n;x++)
	{
		if(x==i) 
		{
			rowFlag = 1;
			continue;
		}
		int colFlag=0;
		for(int y=0;y<n;y++)
		{
			if(y==j)
			{
				 colFlag=1;
				 continue;
			}
			
			// i行以后的行数减 1，j列以后的列数减 1
			result[x-rowFlag][y-colFlag]=arrIni[x][y];
		}
	}
	
	return result;
}

// 伴随矩阵求逆矩阵
double ** inverseMatrix (double ** arr,int n)
{
	// 计算行列式
	double rA = calc2(arr,n);
	
	if(rA == 0) 
	{
		printf("Log[ERROR] - inverseMatrix - 行列式为0，不可逆\n");
		return NULL;
	}
	
	double ** adjointArr = adjoint(arr,n);
	
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			adjointArr[i][j] = adjointArr[i][j]/rA;
		}
	}
	return adjointArr;	
}

// 初等变换求逆矩阵
double ** inverseMatrix2 (double ** arr,int n)
{
	if(rank(arr,n,n)!=n)
	{
		printf("Log[ERROR] - inverseMatrix - 秩小于阶，不可逆\n");
		return NULL;
	}
	
	// 增广矩阵
	double ** big = upDim(malloc(sizeof(double)*n*n*2),n,n*2);

	for(int i = 0;i<n;i++)
	{
		for(int j = 0;j<n;j++)
		{
			big[i][j] = arr[i][j];
		}
		big[i][n+i]=1;
	}

	// 行最简
	bigSimple = simple(big,n,n*2);
	
	double ** result = upDim(malloc(sizeof(double)*n*n),n,n); 
	
	// 单位矩阵初等变换后即逆矩阵
	for(int i = 0;i<n;i++)
	{
		for(int j = 0;j<n;j++)
		{
			result[i][j] = bigSimple[i][n+j];
		}
	}
	recycle(big);
	recycle(bigSimple);
	return result;
}

// 数组副本
double ** copyArr(double ** arr,int row,int col)
{
	int size = sizeof(double)*row*col;
	double ** copyPoint = upDim(malloc(size),row,col);
	memcpy(copyPoint[0],arr[0], size);
	return copyPoint;
}

// 回收数组
void recycle(double ** arr)
{
	free(* arr);
	free(arr);
}
