#include "OpenMP.h"
#include <iostream>
#include <thread>
#include "omp.h"

#define BUFFER_SIZE 1024
#define THREAD_COUNT 2

thread_arg::~thread_arg()
{
    free(data_a);
    free(data_b);
}

thread_arg::thread_arg(coo_data* data_a, coo_data* data_b){
    this->data_a = data_a;
    this->data_b = data_b;
}

coo_data::~coo_data(){
    free(array);
}

coo_data::coo_data(int *data, int size)
{
    this->array = data;
    this->size = size;
}

int* coo_data::get_ptr_at(int i){
    return (i * 3 >= this->size) ? nullptr : (&array[i * 3]);
}

void _merge(thread_arg** args, std::vector<std::vector<int>*>* merged_map_a_ptr, std::vector<std::vector<int>*>* merged_map_b_ptr){
    int j = args[0]->key_size;
    std::vector<std::vector<int>*> merged_map_a = merged_map_a_ptr[0];
    std::vector<std::vector<int>*> merged_map_b = merged_map_b_ptr[0];

    for(int i = 0; i < THREAD_COUNT; i++){
        std::vector<std::vector<int>*> map_a = args[i]->map_a;
        std::vector<std::vector<int>*> map_b = args[i]->map_b;

        // for(int q = 0; q < (*map_a[0]).size(); q++){
        //     std::cout << (*map_a[0])[q] << ", ";
        // }
        // std::cout << std::endl;

        for(int k = 0; k < j; k++){
            while((*map_a[k]).size() != 0){
                merged_map_a[k]->push_back((*map_a[k])[0]);
                merged_map_a[k]->push_back((*map_a[k])[1]);

                (*map_a[k]).erase((*map_a[k]).begin(), (*map_a[k]).begin() + 2);  
            }

            while((*map_b[k]).size() != 0){
                merged_map_b[k]->push_back((*map_b[k])[0]);
                merged_map_b[k]->push_back((*map_b[k])[1]);

                (*map_b[k]).erase((*map_b[k]).begin(), (*map_b[k]).begin() + 2);  
            }
        }
    }
}

void _map_insert(std::vector<int>* map, int x, int y, int v){
    int is_merged = false;

    for(int i = 0; i < map->size(); i += 3){
        
        if ((*map)[i] < x)
            continue;
        
        if ((*map)[i] > x){
            // insert to current index
            map->insert(map->begin() + i, x);
            map->insert(map->begin() + i + 1, y);
            map->insert(map->begin() + i + 2, v);
            return;
        }
        else
        {
            if ((*map)[i + 1] < y)
                continue;

            if ((*map)[i + 1] > y)
            {
                // insert to current index
                map->insert(map->begin() + i, x);
                map->insert(map->begin() + i + 1, y);
                map->insert(map->begin() + i + 2, v);
            }
            else
            {
                // match key found
                (*map)[i + 2] += v;
            }
            return;
        }
    }

    map->push_back(x);
    map->push_back(y);
    map->push_back(v);

}

std::vector<int>* _reduce( std::vector<std::vector<int>*> merged_map_a, std::vector<std::vector<int>*> merged_map_b){
    int j = merged_map_a.size();

    int x, y, v, v_a, v_b;
    int is_merged;

    std::vector<int>* result_map = new std::vector<int>();

    for(int i = 0; i < j; i ++){
        for(int a = 0; a < (*merged_map_a[i]).size() / 2; a++){
            for(int b = 0; b < (*merged_map_b[i]).size() / 2; b++){

                // search in result_map
                x = (*merged_map_a[i])[a * 2];
                y = (*merged_map_b[i])[b * 2];
                
                v_a = (*merged_map_a[i])[a * 2 + 1];
                v_b = (*merged_map_b[i])[b * 2 + 1];

                v = v_a * v_b;
                _map_insert(result_map, x, y, v);

                for(int q = 0; q < (*result_map).size(); q += 3){
                    std::cout << (*result_map)[q] << " " << (*result_map)[q + 1] << " " << (*result_map)[q + 2] << std::endl;
                }
                std::cout << std::endl;
            }
        }
    }

    
    return result_map;
}

void _map(thread_arg *arg)
{
    int j = arg->key_size;

    std::vector<std::vector<int>*>* matrix_map_a;
    matrix_map_a = new std::vector<std::vector<int>*>;
    for(int i = 0; i < j; i++){
        matrix_map_a->push_back(new std::vector<int>());
    }

    std::vector<std::vector<int>*>* matrix_map_b;
    matrix_map_b = new std::vector<std::vector<int>*>;
    for(int i = 0; i < j; i++){
        matrix_map_b->push_back(new std::vector<int>());
    }

    for (int i = 0; i < arg->relative_size; i++)
    {
        // 1. get COO data
        // 2. put to the map

        int *item_a = arg->data_a->get_ptr_at(i);
        int *item_b = arg->data_b->get_ptr_at(i);

        // [i,j,v]
        (*matrix_map_a)[item_a[1]]->push_back(item_a[0]);
        (*matrix_map_a)[item_a[1]]->push_back(item_a[2]);

        // [j,k,v]
        (*matrix_map_b)[item_b[0]]->push_back(item_b[1]);
        (*matrix_map_b)[item_b[0]]->push_back(item_b[2]);
    }

    arg->map_a = *matrix_map_a;
    arg->map_b = *matrix_map_b;
    
}

int main(int argc, char *argv[])
{
    // matrix 1
    // 1 2 0
    // 3 0 5
    // 0 7 0
    // 4 0 0
    // 1 8 2

    // matrix 2
    // 1 4 5 0
    // 0 4 0 1
    // 6 3 9 1

    // matrix_a = i*j
    // matrix_b = j * k
    // key size = j
    int j = 3;
    int key_size = j;

    // mork data
    int coo_array_a[] = {
        0, 0, 1,
        0, 1, 2,
        1, 0, 3,
        1, 2, 5,
        2, 1, 7,
        3, 0, 4,
        4, 0, 1,
        4, 1, 8,
        4, 2, 2};

    int coo_array_b[] = {
        0, 0, 1,
        0, 1, 4,
        0, 2, 5,
        1, 1, 4,
        1, 3, 1,
        2, 0, 6,
        2, 1, 3,
        2, 2, 9,
        2, 3, 1};

    int coo_array_size = *(&coo_array_a + 1) - coo_array_a;
    int processor_count = THREAD_COUNT;
    if (processor_count == 0)
    {
        // not able to get processor_count
        processor_count = 2;
    }

    // omp_set_num_threads(processor_count);

    int size_per_coo_array = (coo_array_size / 3) / processor_count;

    // work devided by number of core
    int x = processor_count;

    int remainder = coo_array_size % processor_count;
    if (remainder != 0)
        x -= 1;

    thread_arg **args = (thread_arg **)malloc(sizeof(thread_arg *) * processor_count);
    for (int i = 0; i < x; i++)
    {
        coo_data *a = new coo_data((coo_array_a + i * size_per_coo_array * 3), size_per_coo_array * 3);
        coo_data *b = new coo_data((coo_array_b + i * size_per_coo_array * 3), size_per_coo_array * 3);
        args[i] = new thread_arg(a, b);
        args[i]->data_a = a;
        args[i]->data_b = b;
        args[i]->relative_size = size_per_coo_array;
        args[i]->key_size = key_size;
    }

    // if the work is not devided evenly.
    if (remainder != 0)
    {
        coo_data *a = new coo_data((coo_array_a + x * size_per_coo_array * 3), size_per_coo_array * 3 + remainder * 3);
        coo_data *b = new coo_data((coo_array_b + x * size_per_coo_array * 3), size_per_coo_array * 3 + remainder * 3);
        args[x] = new thread_arg(a, b);
        args[x]->data_a = a;
        args[x]->data_b = b;
        args[x]->relative_size = size_per_coo_array + remainder;
        args[x]->key_size = key_size;
    }

    // matrix 1 i*j
    // k[0]: [0, 1, 1, 3, 3, 4, 4, 1]
    // ...

    // matrix 2 j*k
    // k[0]: [0, 1, 1, 4, 2, 5]
    // ...

    // MAP

    std::vector<int> *matrix_map_a[j];
    std::vector<int> *matrix_map_b[j];

    omp_set_num_threads(processor_count);
    #pragma omp parallel
    {
    #pragma omp for
        for (int i = 0; i < THREAD_COUNT; i++)
        {
            int thread_idx = omp_get_thread_num();
            args[i]->thread_idx = thread_idx;
            _map(args[thread_idx]);
        }
    }

    #pragma barrier

    std::vector<std::vector<int>*>* merged_map_a_ptr;
    merged_map_a_ptr = new std::vector<std::vector<int>*>;
    for(int i = 0; i < j; i++){
        merged_map_a_ptr->push_back(new std::vector<int>());
    }

    std::vector<std::vector<int>*>* merged_map_b_ptr;
    merged_map_b_ptr = new std::vector<std::vector<int>*>;
    for(int i = 0; i < j; i++){
        merged_map_b_ptr->push_back(new std::vector<int>());
    }

    _merge(args, merged_map_a_ptr, merged_map_b_ptr);
    
    std::vector<int>* result_map = _reduce(*merged_map_a_ptr, *merged_map_b_ptr);
    
    // construct the matrix
    // i * k = 5 * 4 = 20
    int matrix[5][4];
    for(int i = 0; i < 5; i++){
        for(int j = 0; j < 4; j++){
            matrix[i][j] = 0;
        }
    }

    for(int i =0; i < result_map->size(); i+=3){
        matrix[(*result_map)[i]][(*result_map)[i + 1]] = (*result_map)[i + 2];
    }

    for(int i = 0; i < 5; i++){
        for(int j = 0; j < 4; j++){
           std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}
