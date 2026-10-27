#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <cmath> // 新增 cmath 供 exp() 使用
#include <bits/stdc++.h>
#include "saveFitnessHistory.h"
#include <fstream>

using namespace std;

//定義物品結構
/*struct Item{
    string name;
    int weight;
    int value;
    int qty;
};*/

int main(){
    /*vector<Item> base_item = {
        {"A", 1, 6, 10}, {"B", 2, 7, 10}, {"C", 3, 8, 10}, {"D", 4, 9, 10},
        {"E", 5, 10, 10}, {"F", 6, 11, 10}, {"G", 7, 12, 10}, {"H", 8, 13, 10},
        {"I", 9, 14, 10}, {"J", 10, 15, 10}
    };*/
    
    
    // QTS 超參數
    int capacity;
    int num_items;
    int pop_size = 20;
    int generation = 10000;
    double theta = 0.001;
    int num_experiments = 50;
    
    // SA (模擬退火) 超參數
    double T_init = 800.0;       // 初始溫度
    double alpha = 0.9995;     // 降溫係數 (針對 10000 代，0.9995 可以讓 T 穩步降至接近 0)
    
    //將物品攤平
    vector<string> items_names;
    vector<int> items_weights;
    vector<int> items_values;
    
    /*for(const auto& item:base_item){
        for(int i=0; i<item.qty; ++i){
            items_names.push_back(item.name);
            items_weights.push_back(item.weight);
            items_values.push_back(item.value);
        }
    }*/ 
    
    //讀檔、寫入資料 
    ifstream file("200.kp");

    if (!file.is_open()) {
        cout << "Cannot open file!\n";
        return 1;
    }
    
    file >> num_items;
    file >> capacity;

    items_values.resize(num_items);
	items_weights.resize(num_items);
    
    for(int i = 0 ;i<num_items; i++){
    	items_names.push_back(to_string(i + 1));
    	file >> items_values[i] >> items_weights[i];
	}
	
	file.close();
	
	//用DP自己計算最佳解(資料集沒提供最佳解時用) 
    vector<int> dp(capacity+1,0);

	for(int i=0;i<num_items;i++){
	    for(int w=capacity; w>=items_weights[i]; w--){
	        dp[w] = max(dp[w],
	                    dp[w-items_weights[i]]
	                    + items_values[i]);
	    }
	}
	
	cout << "最佳解" <<  dp[capacity] << endl;
    
    random_device rd;
    mt19937 gen(42);
    uniform_real_distribution<double> dist(0.0, 1.0);
    
    //統計變數
    long long sum_best_fitness = 0;   // 記錄 50 次實驗的最佳 Fitness 總和
    int absolute_best_fitness = 0;    // 記錄 50 次實驗中的最佳 Fitness
    vector<int> absolute_best_solution;
    int absolute_best_weight = 0; 
    int best_exp = 0; //紀錄50次實驗中最佳解出現在哪一次 
    
    // 建立二維陣列來存放所有的歷史紀錄 
    // 大小為：num_experiments 乘以 generation的，初始值為 0.0
    vector<vector<double>> allHistory(num_experiments, vector<double>(generation, 0.0));
    
    for(int exp = 0; exp < num_experiments; ++exp) { 
    
        //演算法變數
        vector<double> item_prob(num_items, 0.5);
        vector<int> global_best_solution(num_items, 0);
        int global_best_fitness = 0;
        int global_best_weight = 0;
        
        // 每次實驗重置溫度
        double T = T_init; 
        
        //演算法主迴圈
        for(int iter=0; iter<generation; ++iter){
            vector<vector<int>> population(pop_size, vector<int>(num_items, 0));
            vector<int> fitnesses(pop_size, 0);
            vector<int> weights(pop_size, 0);
            
            for(int i = 0; i<pop_size; ++i){
                int total_weight = 0;
                int total_value = 0;
                
                // 觀測並計算 fitness
                for(int j = 0; j < num_items; ++j){
                    if(dist(gen) < item_prob[j]){
                        population[i][j] = 1;
                        total_value += items_values[j];
                        total_weight += items_weights[j];
                    }
                }
                
                weights[i] = total_weight;
                //懲罰機制
                fitnesses[i] = (total_weight > capacity) ? 0 : total_value;
            }
        
            // 找出當代最佳與最差解
            int best_idx = 0, worst_idx = 0;
            for(int i = 1; i < pop_size; ++i){
                if(fitnesses[i] > fitnesses[best_idx]) best_idx = i;
                if(fitnesses[i] < fitnesses[worst_idx]) worst_idx = i;
            }
            
            //計算delta
			int delta  =  global_best_fitness - fitnesses[best_idx];
        
            //更新歷史最佳解
            if(delta < 0){
                global_best_fitness = fitnesses[best_idx];
                global_best_solution = population[best_idx];
                global_best_weight = weights[best_idx];
            } else {
            	//計算接受差解的機率 
            	double prob = std::exp(-(double)delta / T);
            	if (dist(gen) < prob){
            		global_best_fitness = fitnesses[best_idx];
                	global_best_solution = population[best_idx];
                	global_best_weight = weights[best_idx];
				}
			}
            
            //將每代的global_best_fitness存放至allHistory
            allHistory[exp][iter] = global_best_fitness;
            
            
            // 更新量子機率
            for(int j = 0; j < num_items; ++j){
	            int b_bit = global_best_solution[j];
	            int w_bit = population[worst_idx][j];
	            
	            if (b_bit == 1 && w_bit == 0){
	                item_prob[j] = min(0.95, item_prob[j] + theta);
	            } else if (b_bit == 0 && w_bit == 1){
	                item_prob[j] = max(0.05, item_prob[j] - theta);
	            }
	        }
            
            // 每代結束後進行降溫
            T *= alpha;
            
        } //世代迴圈結束 
        cout << "Exp " << exp + 1 << " Best Fitness: " << global_best_fitness << endl;
        sum_best_fitness += global_best_fitness;
        
        // 記錄50次實驗中的整體歷史最佳解
        if(global_best_fitness > absolute_best_fitness){
            absolute_best_fitness = global_best_fitness;
            absolute_best_solution = global_best_solution;
            absolute_best_weight = global_best_weight;
            best_exp = exp + 1;
        }
    } 
    
    //輸出結果與統計
    double avg_best_fitness = (double)sum_best_fitness / num_experiments;
    cout << "\nCompleted " << num_experiments << " Experiments" << endl;
    cout << "Average Best Fitness: " << avg_best_fitness << endl;
    cout << "\nAbsolute Best Record" << endl;
    cout << "Best experiment: " << best_exp << endl;
    cout << "Best value: " << absolute_best_fitness << endl;
    cout << "Weight: " << absolute_best_weight << " / " << capacity << endl; 
    
    /*cout << "\n物品選擇數量統計:" << endl;
    int current_idx = 0;
    for(const auto& item:base_item){
        int selected_count = 0;
        for(int i = 0; i < item.qty; ++i){
            selected_count += absolute_best_solution[current_idx++];
        }
        cout << item.name << ": " << selected_count << "個" << endl;
    }*/
    
    saveFitnessHistory("SAQTS_FitnessHistory_200_T1000-3.csv", allHistory);
    
    return 0;
}
