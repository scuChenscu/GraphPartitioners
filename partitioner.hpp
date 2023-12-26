//
// Created by muzongshen on 2021/9/23.
//

#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

#include "util.hpp"

using namespace std;
/**
 * Partitioner，分区器，其他所有分区算法都继承该类
 */
class Partitioner
{
protected:
    Timer total_time;
    // 输出文件流对象
    ofstream edge_fout;
    ofstream node_fout;
    // 负载不均衡
    // 边割率
    double edge_cut;
    double load_balance;
public:
    // 分区方法，纯虚函数，派生类必须重写该虚函数
    virtual void split() = 0;
    void set_write_files(const string &basefilename,const string &method, int pnum){
        // TODO 传入文件名，这里需要修改成对应的文件名
        edge_fout.open(partitioned_name(basefilename,method,pnum));
        node_fout.open(nodesave_name(basefilename,method,pnum));
    }
    // 保存每个vertex及其对应的分区， v_id p_id
    void save_vertex(vid_t v, int p_id)
    {
        node_fout << v << " "<< p_id <<endl;
    }
    // 保存每条edge分区，from to p_id, 应该是使用了点分割和边分割的区别
    void save_edge(vid_t from, vid_t to, int p_id)
    {
        edge_fout<<from<<" "<<to<<" "<<p_id<<endl;
    }

    double calculate_edge_cut()
    {
        edge_cut = 0;
        for (int i = 0; i < p; i++)
        {
            edge_cut += (double)subg_vids[i].size() / (double)num_vertices;
        }
        return edge_cut;
    }

    double calculate_load_balance()
    {
        load_balance = 0;
        for (int i = 0; i < p; i++)
        {
            load_balance += (double)subg_vids[i].size() / (double)num_vertices;
        }
        return load_balance;
    }
};
