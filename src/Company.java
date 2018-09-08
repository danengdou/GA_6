import java.io.*;
import java.util.*;

public class Company {

	private static int rows;				//数据的行数
	private static int column;			//数据的列数
	private static String dataPath;		//输入路径
	private static String resultPath;	//输出路径
	private static int user_resource;    //用户自定义资源数量
	private static int popSize;          //种群数量
	private static int maxgens;          //迭代次数
	private static double pxover;        //交叉概率
	private static double pmultation;    //变异概率
	private static int g_max;            //优化重合数的迭代次数
	private static int pen_more;         //超过资源数的惩罚
	private static int pen_re;           //实施和计划重合的惩罚
	private static int avg_re_num;       //均衡资源的数量
	private static int[][] realArray;    //数据读取
	private static int key;              //计划提前与计划平均相冲突，此处为切换开关，key=1代表提前，key=0代表平均
	private static double pmultationAdv; //计划提前程度,当key=1时，设置该值，取值范围为0.5~1，越大越提前，越小越平均
	private static int[] lockRows;       //锁定行
	private static int[] lockColumns;    //锁定列
	private static ArrayList lockRowsList;       //锁定行集合
	private static ArrayList lockColumnsList;    //锁定列集合
	
	
	public Company() {
		this.rows = 0;
		this.column = 0;
		//this.dataPath = "E:\\Java\\workspace\\GA5\\data.csv";
		//this.resultPath = "E:\\Java\\workspace\\GA5\\result.csv";
		this.dataPath = "C:\\Users\\hantao5\\Desktop\\GA_6\\Book1.csv";
		this.resultPath = "C:\\Users\\hantao5\\Desktop\\GA_6\\result.csv";

		this.user_resource = 6;
		this.popSize = 50;
		this.maxgens = 10000;
		this.pxover = 0.8;
		this.pmultation = 0.2;
		this.g_max = 500;
		this.pen_more = 10000;
		this.pen_re = 1;
		this.key = 0;
		this.pmultationAdv = 1;
		this.lockRows = new int[]{};
		this.lockColumns = new int[]{};
	}
	 
	public static void main(String[] args) {
		Company company = new Company();
		company.process();
	}

	//主要逻辑
	private void process() {
		ArrayList ini_List = readData();         		//简化形式读取excel中的数据-计划方案
		realArray = readData2();         		        //完整形式读取excel中的数据-计划方案
		int[] firstColumn = saveFirstColumn(realArray); //保存第1列的值 -- 解决第1列值有0报异常
		System.out.println("----------计划方案为---------");
		System.out.println(ini_List);					//打印计划Arraylist
		setRealArrayIsNotZero(realArray);//数组中第一列的值暂时置为1
		int len_of_gene = ini_List.size();				//染色体的长度
		int numAllRes = addAllRes(realArray);           //计算总资源数量
		int[] plan  = changeArrToGroup(ini_List);		//将Arraylist转化为数组
		avg_re_num = avgReNum(numAllRes);			    //均衡资源数量
		int count = 0;
		lockRowsList = arrayToList(lockRows);           //将锁定行数组转为集合
		lockColumnsList = arrayToList(lockColumns);     //将锁定列数组转为集合
		int [][]farm = createFarmAvg(popSize,plan);                          //0、生成初始解
		while(count < maxgens){
			int[][] paixu_farm = choPaiXu(farm, plan);                       //1、染色体排序
			choXuanZe(paixu_farm, plan);									 //2、染色体选择
			choCross(paixu_farm, plan, len_of_gene);           				 //3、染色体交叉
			if (key == 1) {
				choMutateAdv(paixu_farm, plan, len_of_gene);                 //4、染色体变异 - 计划提前好
			} else {
				choMutateAvg(paixu_farm, plan, len_of_gene);                 //4、染色体变异 - 计划平均好
			}
			farm = paixu_farm;
			count=count+1;
		}
		printResult(farm, plan, firstColumn);
	}

	//数组转集合
	private ArrayList arrayToList(int[] array) {
		ArrayList<Integer> a = new ArrayList<Integer>();
		if (array.length > 0) {
			for (int i=0;i<array.length;i++) {
				a.add(array[i]);
			}
		}
		return a;
	}
	//计算总资源个数
	private int addAllRes(int[][] realArray) {
		int count = 0;
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < column; j++) {
				count += realArray[i][j];
			}
		}
		return count;
	}

	public int[] saveFirstColumn(int[][] realArray) {
		int[] firstColumn = new int[rows];
		for (int i=0; i<rows; i++) {
			firstColumn[i] = realArray[i][0];
		}
		return firstColumn;
	}

	private void setRealArrayIsNotZero(int[][] realArray) {
		for (int i=0; i<rows; i++) {
			realArray[i][0] = 1;
		}
	}
	private void printResult(int[][] farm, int[] plan, int[] firstColumn) {
		int[][] paixu_farm = choPaiXu(farm, plan);                    
		int[] best_result = paixu_farm[0];//存储最优路径
		int[] farm_fitness2 = calFitness(paixu_farm,plan);//计算染色体的适应度
		int[] compareData = groupData(best_result);//将完整数据拼成一维数组
		int[][] result = createData(best_result, compareData, firstColumn);	//根据一维数组还原矩阵
		 if(farm_fitness2[0] < pen_more){
			 FileWriter fw = null;
			 BufferedWriter bw = null;
			 System.out.println("--------------------------------------------------------------------");
			 System.out.println("--------满足条件的解-------");
			 System.out.println(Arrays.toString(best_result));
			 System.out.println("--------满足条件的解矩阵形式-------");
			 try {
				 fw =new FileWriter(resultPath);
				 bw=new BufferedWriter(fw);
				 for(int i=0;i<result.length;i++){
					 bw.write(Arrays.toString(result[i]).replace("[", "").replace("]", "")+"\r\n");
				 }
				 System.out.println(Arrays.deepToString(result));
			} catch (IOException e) {
				e.printStackTrace();
			}
			 finally
		        {
		            try {
						bw.close();
						fw.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
		            
		        }
		 } else {
			 //输出不可行解
			 FileWriter fw = null;
			 BufferedWriter bw = null;
			 System.out.println("----------------------------");
			 System.out.println("可能不存在满足条件的解，最优的非可行解为:");
			 System.out.println(Arrays.toString(best_result));
			 System.out.println("-----资源数量超过阈值的时间段是----");
			 
			 printNoRes(best_result);                                    //输出不可行解的资源数和时间
			 
			 System.out.println("----非可行解的矩阵形式----");
			 try {
				 fw =new FileWriter(resultPath);
				 bw=new BufferedWriter(fw);
				 for(int i=0;i<result.length;i++){
					 bw.write(Arrays.toString(result[i]).replace("[", "").replace("]", "")+"\r\n");
				 }
				 System.out.println(Arrays.deepToString(result));
			} catch (IOException e) {
				e.printStackTrace();
			}
			 finally
		        {
		            try {
						bw.close();
						fw.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
		            
		        } 
		 }
	}

	//根据一维数组还原矩阵
	private int[][] createData(int[] best_result, int[] compareData, int[] firstColumn) {
		 int[][] result = new int[rows][column];
		 int enter = 0;
		 for(int i=0;i<best_result.length;i++){
			 if(best_result[i]==0){
					enter = enter + 1;
					result[enter-1][best_result[i]] = compareData[i];
				}else{
					result[enter-1][best_result[i]] = compareData[i];
				}
		 }
		 for (int i=0;i<rows;i++) {
			 result[i][0] = firstColumn[i];
		 }
		return result;
	}

	//将完整数据拼凑成一维数组
	private int[] groupData(int[] best_result) {
		int[] array = new int[best_result.length];
		int k = 0;
		for (int i = 0; i < realArray.length; i++) {
			for (int j = 0; j < realArray[0].length; j++) {
				if (realArray[i][j] != 0) {
					array[k++] = realArray[i][j];
				}
			}
		}
		return array;
	}

	//输出超过阈值的资源数及时间
	private void printNoRes(int[] best_result) {
		 int hang = -1;
		 int[] everyResource = new int[column];
		 for (int i = 0; i < best_result.length; i++) {
			if (best_result[i] == 0) {
				hang ++;
			}
			if (best_result[i] != 0) {
				int lie = best_result[i];
				while ((lie < column) && (realArray[hang][lie] == 0)) {
					++ lie;
				}
				everyResource[best_result[i]] += realArray[hang][lie];
			}
	 	 }
		 for (int i = 1; i < everyResource.length; i++) {
			if (everyResource[i] > user_resource) {
				System.out.println("第"+i+"个时间段资源数量为"+everyResource[i]+"超过了阈值"+user_resource);
			}
		 }
	}

	//简化形式读取数据
	private static ArrayList readData() {
		 File file = new File(dataPath);
		 BufferedReader br = null;
	     FileReader fb = null;
	     ArrayList ini_List = new ArrayList(); 
	     String[] temp=null;;
	        try
	        {
	            fb = new FileReader(file);
	            br = new BufferedReader(fb);
	            String str = null;
	            while ((str = br.readLine()) != null)
	            {
	            	 rows++;
		             ini_List.add(0);
	            	 temp = str.split(","); 
	            	 for(int j=1;j<temp.length;j++){
	            		 //分割字符串
	 	                if ( Integer.parseInt( temp[j] )!=0 ){
	 	                	ini_List.add(j);
	 	                }  
	 	            }  
	            }
	        }
	        catch (IOException e)
	        {
	            e.printStackTrace();
	        }
	        finally
	        {
	            try {
					br.close();
					fb.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
	            
	        }
	        column = temp.length;   //列数，即时间段数
		 return ini_List;
	}

	//完整形式读取数据
	private static int[][] readData2() {
		 File file = new File(dataPath);
		 BufferedReader br = null;
	     FileReader fb = null;
	     int[][] realArray = new int[rows][column];
	     String[] temp=null;;
	        try
	        {
	            fb = new FileReader(file);
	            br = new BufferedReader(fb);
	            String str = null;
	            for (int i = 0; i < rows; i++) {
	            	str = br.readLine();
	            	temp = str.split(","); 
					for (int j = 0; j < column; j++) {
						realArray[i][j] = Integer.parseInt( temp[j] );
					}
				}
	        }
	        catch (IOException e)
	        {
	            e.printStackTrace();
	        }
	        finally
	        {
	            try {
					br.close();
					fb.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
	            
	        }
		 return realArray;
	}
	
	//备份染色体
	private static int[][] copy_farm(int len_of_gene,int[][] paixu_farm){
		int[][] paixu_farm2 = new int[popSize][len_of_gene];
		for (int i=0;i<popSize;i++){
			for(int j=0;j<len_of_gene;j++){
				paixu_farm2[i][j] = paixu_farm[i][j];
			}
		}
		return paixu_farm2;
	}
	
    //把ArrayList转成数组
  	private static int[] changeArrToGroup(ArrayList order){
  		int width = order.size();//染色体的长度
  		int[] save_list = new int[width];
  			   for (int j = 0; j <width; j++) {
  				  //
			       save_list[j]= Integer.parseInt(order.get(j).toString());
  			   }
  		return save_list;
      }
  	
  	//计算均衡资源数量（（总资源数-0的个数）/（时间段数-1））
	private static int avgReNum(int len_of_gene) {
		int tempnum = (int) Math.round( ((len_of_gene - rows) * 1.0) / (column - 1));
		return tempnum;
	}
	
	//计算均衡资源数的适应度
	private static int avgFit(int[] chor){
		int[] everyResource = new int[column];
		int hang = -1;
		for (int j = 0; j < chor.length; j++) {
			if (chor[j] == 0) {
				hang ++;
			}
			if (chor[j] != 0) {
				int lie = chor[j];
				while ((lie < column) && (realArray[hang][lie] == 0)) {
					++ lie;
				}
				everyResource[chor[j]] += realArray[hang][lie];
			}
		}
		int temp = 0;
		for(int i = 1; i < everyResource.length; i++) {
			temp += Math.abs(everyResource[i] - avg_re_num);
		}
		return temp;
	}
	
	//计算资源数超过阈值的适应度
	private static int resFit(int[] chor){
		int[] everyResource2 = new int[column];
		int hang = -1;
		for (int j = 0; j < chor.length; j++) {
			if (chor[j] == 0) {
				hang ++;
			}
			if (chor[j] != 0) {
				int lie = chor[j];
				while ((lie < column) && (realArray[hang][lie] == 0)) {
					++ lie;
				}
				everyResource2[chor[j]] += realArray[hang][lie];
			}
		}
		int temp = 0;
		for(int i = 1; i < everyResource2.length; i++) {
			if (everyResource2[i] > user_resource) {
				temp += everyResource2[i] - user_resource;
			}
		}
		return temp*pen_more;
	}
	
	//计算染色体适应度
  	private static int[] calFitness(int[][] farm ,int[] plan){
  		int num_farm = farm.length;
  		int []farm_fitness = new int[num_farm];
  		int fitness = 0;
  		for (int i=0;i<num_farm;i++){
  			fitness += resFit(farm[i]) + avgFit(farm[i]);
  			farm_fitness[i] = fitness;
  		}
		return farm_fitness;
  	}
   
  	//染色体排序
  	private static int[][] choPaiXu(int[][] farm, int[] plan){
		 int[] farm_fitness = calFitness(farm,plan);//计算染色体的适应度
		 int[] farm_fitness2= Arrays.copyOf(farm_fitness, farm_fitness.length);//复制染色体的适应度数组
		 Arrays.sort(farm_fitness2);//对适应值从小到大排序
		 int[] move_refit = remove_refitness(farm_fitness2);//去掉适应度中重复的值
		 int[] fit_index = index_paixu(farm_fitness,move_refit);//查找排序后适应度所对应的染色体的下标，并将染色体按适应度排序
		 int[][] paixu_farm = farm_paixu(plan.length,farm,fit_index);//将染色体按适应度排序
		 return paixu_farm;
  	}
   
  	//染色体选择
  	private static void choXuanZe(int[][] paixu_farm, int[] plan) {
  		int k = 2;//用最优的染色体替换最差的染色体的条数
  		 for (int i=0;i<k;i++){
  				for(int j=0;j<plan.length;j++){
  					paixu_farm[popSize-i-1][j] = paixu_farm[i][j];
  				}
  		 }
  	}
 
  	//染色体交叉
    private static void choCross(int[][] paixu_farm, int[] plan, int len_of_gene){
		int[][] paixu_farm2 = copy_farm(len_of_gene, paixu_farm);            //备份排序后染色体
    	//3、-----------------交叉
		 //---------根据染色体的长度确定交叉检验的位数
		 int cross_num = 1;
		 if (len_of_gene<10){
			  cross_num = 2;
		 }else if (len_of_gene>10 && len_of_gene<20){
			  cross_num = 3;
		 }else if(len_of_gene>20 && len_of_gene<40){
			  cross_num = 4;
		 }else if(len_of_gene>40 && len_of_gene<80){
			  cross_num = 6;
		 }else{
			  cross_num = 10;
		 }
		//---------计算总的交叉次数
		 int sum_cross_num = (int) (pxover * popSize * cross_num ) ;
	    //---------交叉
		 int A1,A2,cross_position,temp1,temp2;
		 for (int i=0;i<sum_cross_num;i++){
			 Random rand = new Random();
			  A1 = rand.nextInt(popSize-1);
			  A2 = rand.nextInt(popSize-1);
			  cross_position = rand.nextInt(len_of_gene-1);
			  //保证cross_position不在锁定范围
			  if (paixu_farm[A1][cross_position] == 0 || paixu_farm[A2][cross_position] == 0) {
				 continue; 
			  } else {
				  if (lockRowsList.contains(whatRow(cross_position,plan)) 
						  || lockColumnsList.contains(paixu_farm[A1][cross_position]+1)
						  || lockColumnsList.contains(paixu_farm[A2][cross_position]+1) ) {
					  continue; 
				  }
			  }
			  temp1 = paixu_farm[A1][cross_position];
			  temp2 = paixu_farm[A2][cross_position];
			  paixu_farm[A1][cross_position] = temp2;
			  paixu_farm[A2][cross_position] = temp1;
		 }
		 int[] farm_cross_fitness = calFitness(paixu_farm, plan);//计算交叉后染色体的适应度
		 int[] farm_fitness2 = calFitness(paixu_farm2, plan);//计算交叉后前染色体的适应度
		 //----------交叉检验(若交叉后染色体适应度比原先的差，则不采用新染色体)
		 for (int i=0;i<popSize;i++){
			 if(farm_fitness2[i]<farm_cross_fitness[i]){
				 paixu_farm[i] = paixu_farm2[i];
			 }
		 }

    }
	//计算一个数位于第几行
    private static int whatRow(int num, int[] plan) {
    	int k = 0;
    	for (int i=0;i<num;i++) {
    		if (plan[i] == 0) {
    			k++;
    		}
    	}
    	return k;
    }
    //染色体变异
    private static void choMutateAdv(int[][] paixu_farm, int[] plan, int len_of_gene){
    	 int[][] paixu_farm3 = copy_farm(plan.length, paixu_farm);            //备份交叉后染色体
    	 //---------根据染色体的长度确定交叉检验的位数
		 int mutate_num = 1;
		 if (len_of_gene<10){
			 mutate_num = 2;
		 }else if (len_of_gene>=10 && len_of_gene<20){
			 mutate_num = 3;
		 }else if(len_of_gene>=20 && len_of_gene<40){
			 mutate_num = 4;
		 }else if(len_of_gene>=40 && len_of_gene<60){
			 mutate_num = 6;
		 }else if(len_of_gene>60 && len_of_gene<80){
			 mutate_num = 8;
		 }else{
			 mutate_num = 10;
		 }
		//---------计算总的变异次数
		 int sum_mutate_num = (int) (pmultation * popSize * mutate_num ) ;
	    //---------变异
		 int B1,mutate_position;
		 for (int i=0;i<sum_mutate_num;i++){
			 Random rand = new Random();
			 B1 = rand.nextInt(popSize-1);
			 mutate_position = rand.nextInt(len_of_gene-1);
			//保证mutate_position不在锁定范围
			if (lockRowsList.contains(whatRow(mutate_position,plan))|| lockColumnsList.contains(paixu_farm[B1][mutate_position]+1)) {
			     continue; 
			}
			  if(plan[mutate_position] == 0){
				 paixu_farm[B1][mutate_position]=0;
				}else{
					double random2;
					random2 = Math.random();
					double point = plan[mutate_position-1] +  (plan[mutate_position]-plan[mutate_position-1])/2.0;
					if (point - 0.5 > plan[mutate_position-1]) {
						int turePoint = (int)Math.floor(point);//实际分割值
						if (random2 <=pmultationAdv){
							int temp1 = rand.nextInt(turePoint-plan[mutate_position-1])+1+plan[mutate_position-1];
							if (lockColumnsList.contains(temp1+1)) {
								continue;
							} else {
								paixu_farm[B1][mutate_position] = temp1;
							}
							
						}else{
							int temp2 = rand.nextInt(plan[mutate_position]-turePoint)+1+turePoint;
							if (lockColumnsList.contains(temp2+1)) {
								continue;
							} else {
								paixu_farm[B1][mutate_position] = temp2;
							}
						}
					} else {
						paixu_farm[B1][mutate_position] = plan[mutate_position];
					}
				}
		 }

		 int[] farm_mutate_fitness = calFitness(paixu_farm, plan);//计算变异后染色体的适应度
		 int[] farm_fitness = calFitness(paixu_farm3, plan);//计算变异前染色体的适应度

		 //----------变异检验(若变异后染色体适应度比原先的差，则不采用新染色体)
		 for (int i=0;i<popSize;i++){
			 if(farm_fitness[i]<farm_mutate_fitness[i]){
				 paixu_farm[i] = paixu_farm3[i];
			 }
		 }
    }
    
  	//去掉适应度中重复的值
  	private static int[] remove_refitness(int[] farm_fitness2){
	  	Set<Integer> set = new TreeSet<Integer>();//新建一个set集合保存不重复的适应值
		for (int i : farm_fitness2) {
			set.add(i);
		}
	    Integer[] farm_fitness3 = set.toArray(new Integer[0]);// 数组的包装类型不能转 只能自己转；把Integer转为为int数组；
		int[] farm_fitness4 = new int[farm_fitness3.length];
		for (int i = 0; i < farm_fitness4.length; i++) {
			farm_fitness4[i] = farm_fitness3[i];
		}
		return farm_fitness4;
  	}

  //获得染色体排序的下标
  	private static int[] index_paixu(int[] farm_fitness,int[] move_refit){  
	  	int[] fit_index = new int[farm_fitness.length];
		 int k = 0;
		for(int i=0;i<move_refit.length;i++){
			for(int j=0;j<farm_fitness.length;j++){
				if(move_refit[i]==farm_fitness[j]){
					fit_index[k] = j;
					k = k+1;
				}
			}
		}
		return fit_index;
  	}
  	
  //将染色体按适应度排序
  	private static int[][] farm_paixu(int len_of_gene,int[][] farm,int[] fit_index){ 
	  	int[][] paixu_farm = new int[popSize][len_of_gene];
		 for (int i=0;i<popSize;i++){
			 for (int j=0;j<len_of_gene;j++){
				 paixu_farm[i][j] = farm[fit_index[i]][j];
			 }
		 }
		 return paixu_farm;
	 }

	//生成初始种群
	private static int[][] createFarmAvg(int popSize,int[] plan){
		Random rand = new Random();
		int [][] farm = new int[popSize][plan.length];
		for (int i=0;i<popSize;i++){
			int h = 0;
			for (int j=0;j<plan.length;j++){
				if (lockRowsList.contains(h) || lockColumnsList.contains(plan[j]+1)) {
					farm[i][j]=plan[j];
					if(plan[j] == 0) {
						h++;
					}
				} else {
					if(plan[j] == 0){
						farm[i][j]=0;
						h++;
					}else{
						int temp = rand.nextInt(plan[j]-plan[j-1]) + 1 + plan[j-1];
						if (lockColumnsList.contains(temp+1)) {
							farm[i][j] = plan[j];
						} else {
							farm[i][j] = temp;
						}
					}
				}
			}
		}
		return farm;
	}
	
    //染色体变异 - 计划平均好
    private static void choMutateAvg(int[][] paixu_farm, int[] plan, int len_of_gene){
    	 int[][] paixu_farm3 = copy_farm(plan.length, paixu_farm); 
		 int mutate_num = 1;
		 if (len_of_gene<10){
			 mutate_num = 2;
		 }else if (len_of_gene>=10 && len_of_gene<20){
			 mutate_num = 3;
		 }else if(len_of_gene>=20 && len_of_gene<40){
			 mutate_num = 4;
		 }else if(len_of_gene>=40 && len_of_gene<60){
			 mutate_num = 6;
		 }else if(len_of_gene>60 && len_of_gene<80){
			 mutate_num = 8;
		 }else{
			 mutate_num = 10;
		 }
		 int sum_mutate_num = (int) (pmultation * popSize * mutate_num ) ;
		 int B1,mutate_position;
		 for (int i=0;i<sum_mutate_num;i++){
			 Random rand = new Random();
			 B1 = rand.nextInt(popSize-1);
			 mutate_position = rand.nextInt(len_of_gene-1);
			//保证mutate_position不在锁定范围
			  if (lockRowsList.contains(whatRow(mutate_position,plan))|| lockColumnsList.contains(paixu_farm[B1][mutate_position]+1)) {
			 	     continue; 
			  }
			  if(plan[mutate_position] == 0){
				 paixu_farm[B1][mutate_position]=0;
				}else{
					int temp3 = rand.nextInt(plan[mutate_position]-plan[mutate_position-1]) + 1 + plan[mutate_position-1];
					if (lockColumnsList.contains(temp3+1)){
						continue;
					} else {
						paixu_farm[B1][mutate_position] = temp3;
					}
				}
		 }
		 int[] farm_mutate_fitness = calFitness(paixu_farm, plan);
		 int[] farm_fitness = calFitness(paixu_farm3, plan);
		 for (int i=0;i<popSize;i++){
			 if(farm_fitness[i]<farm_mutate_fitness[i]){
				 paixu_farm[i] = paixu_farm3[i];
			 }
		 }
    }
}
