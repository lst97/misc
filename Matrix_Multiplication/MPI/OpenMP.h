// Using MapReduce - code for node

#include <vector>
// matrix_a, matrix_b will be reduce to COOridinate format from server
// store indices of the matrix that have non zero values

// data formate need to be data_a: [[i, j, v]...[i, j, v]]
// data formate need to be data_b: [[j, k, v]...[j, k, v]]
struct coo_data{
    private:
    int* array;

    public:
    coo_data(int *data, int size);
    // int* operator[](int i);
    int* get_ptr_at(int i);
    int size;

    ~coo_data();
};

struct thread_arg
{
    coo_data* data_a;
    coo_data* data_b;
    std::vector<std::vector<int>*> map_a;
    std::vector<std::vector<int>*> map_b;
    int thread_idx;

    int relative_size;
    int key_size;

    public:
    thread_arg(coo_data* data_a, coo_data* data_b);

    ~thread_arg();
};