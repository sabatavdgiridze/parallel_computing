#include <cuda_runtime.h>
#include <cuda.h>



#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

#define TRUE 1
#define FALSE 0

void swap_int(int * first, int * second) {
    int temp = * first;

    * first = * second;
    * second = temp;
}

struct node {
    int key;
    node * left;
    node * right;
};


struct node * create_tree(int * array, int length) {
    
    struct node * root = (struct node *) malloc(sizeof(struct node));
    
    root->left = NULL;
    root->right = NULL;
    root->key = array[0];

    int idx;
    for(idx = 1; idx < length; idx++) {
        struct node * current = root;
        while (TRUE) {
            if (current->key < array[idx]) {
                if (current->right == NULL) {
                    current->right = (struct node *) malloc(sizeof(struct node));

                    current->right->key = array[idx];

                    current->right->left = NULL;
                    current->right->right = NULL;

                    break;
                } else {
                    current = current->right;
                }
            } else {
                if (current->left == NULL) {
                    current->left = (struct node *) malloc(sizeof(struct node));

                    current->left->key = array[idx];

                    current->left->left = NULL;
                    current->left->right = NULL;

                    break;
                } else {
                    current = current->left;
                }
            }
        }
    }
    return root;
}

int flatten_tree_rec(int position, struct node * root, int * keys, int * first_child, int * second_child);

void flatten_tree(struct node * root, int * keys, int * first_child, int * second_child) {
    flatten_tree_rec(0, root, keys, first_child, second_child);
}

int flatten_tree_rec(int position, struct node * root, int * keys, int * first_child, int * second_child) {
    keys[position] = root->key;
    

    int new_position = position + 1;

    if (root->left != NULL) {
        first_child[position] = new_position;
        new_position = flatten_tree_rec(new_position, root->left, keys, first_child, second_child);
    } else {
        first_child[position] = -1;
    }

    if (root->right != NULL) {
        second_child[position] = new_position;
        new_position = flatten_tree_rec(new_position, root->right, keys, first_child, second_child);
    } else {
        second_child[position] = -1;
    }
    return new_position;
}


int tree_search(int key, struct node * root) {
    if (root == NULL) {
        return FALSE;
    } else {
        if (root->key == key) {
            return TRUE;
        } else if (root->key < key) {
            return tree_search(key, root->right);
        } else {
            return tree_search(key, root->left);
        }
    }
}

int flatten_tree_search(int key, int * keys, int * first_child, int * second_child) {
    int root = 0;
    while(TRUE) {
        if (keys[root] == key) {
            return TRUE;
        } else if (keys[root] < key) {
            if (second_child[root] == -1) {
                break;
            } else {
                root = second_child[root];
            }
        } else {
            if (first_child[root] == -1) {
                break;
            } else {
                root = first_child[root];
            }
        }
    }
    return FALSE;
}


__global__ void cuda_flatten_tree_search(int num_of_queries, int * search_keys, int * keys, int * first_child, int * second_child, int * answers) {
    
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx < num_of_queries) {
        int root = 0;
        while(TRUE) {
            if (keys[root] == search_keys[idx]) {
                answers[idx] = TRUE;
                break;
            } else if (keys[root] < search_keys[idx]) {
                if (second_child[root] == -1) {
                    answers[idx] = FALSE;
                    break;
                } else {
                    root = second_child[root];
                }
            } else {
                if (first_child[root] == -1) {
                    answers[idx] = FALSE;
                    break;
                } else {
                    root = first_child[root];
                }
            }
        }
    }
} 

void print_array(int * array, int length) {
    int idx;
    for (idx = 0; idx < length; idx++) {
        printf("%4d ", array[idx]);
    }
    printf("\n");
}

void shuffle(int * array, int length) {
    int idx;
    for (idx = length - 1; idx >= 1; idx--) {
        int j_idx = rand() % (idx + 1);
        swap_int(&array[idx], &array[j_idx]);
    }
}

void allocate_device(void ** handle, int length) {
    if (cudaMalloc(handle, sizeof(int) * length) != cudaSuccess) {
        exit(EXIT_FAILURE);
    }
}

void copy_to_device(int * device_handle, int * host_handle, int length) {
    if (cudaMemcpy(device_handle, host_handle, sizeof(int) * length, cudaMemcpyHostToDevice) != cudaSuccess) {
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char ** argv)
{
    int N = atoi(argv[1]);
    int Seed = atoi(argv[2]);

    srand(Seed);

    int * array = (int *) malloc(2 * N * sizeof(int));
    
    int idx;
    for (idx = 0; idx < 2 * N; idx++) {
        array[idx] = idx;
    }

    shuffle(array, 2 * N);
    if (N <= 20) print_array(array, N);

    int * keys = (int *) malloc(N * sizeof(int));
    int * first_child = (int *) malloc(N * sizeof(int));
    int * second_child = (int *) malloc(N * sizeof(int));

    struct node * root = create_tree(array, N);

    shuffle(array, 2 * N);
    if (N <= 20) print_array(array, N);


    flatten_tree(root, keys, first_child, second_child);



    int * answers_tree = (int *) malloc(N * sizeof(int));
    int * answers_flattened_tree = (int *) malloc(N * sizeof(int));
    
    clock_t start;
    clock_t finish;

    start = clock();
    for (idx = 0; idx < N; idx ++) {
        answers_tree[idx] = tree_search(array[idx], root);
    }
    finish = clock();

    printf("searching in tree data structure took %d miliseconds\n", (finish - start) * 1000 / CLOCKS_PER_SEC);

    start = clock();
    for (idx = 0; idx < N; idx ++) {
        answers_flattened_tree[idx] = flatten_tree_search(array[idx], keys, first_child, second_child);
    }
    finish = clock();

    printf("searching in flattened tree data structure took %d miliseconds\n", (finish - start) * 1000 / CLOCKS_PER_SEC);

    if (N <= 20) print_array(answers_tree, N);
    if (N <= 20) print_array(answers_flattened_tree, N);

    int * answers_flattened_tree_cuda = (int *) malloc(N * sizeof(int));

    int * device_keys;
    int * device_search_keys;
    int * device_first_child;
    int * device_second_child;
    int * device_answers;

    allocate_device((void**)&device_keys, N);
    allocate_device((void**)&device_search_keys, N);
    allocate_device((void**)&device_first_child, N);
    allocate_device((void**)&device_second_child, N);
    allocate_device((void**)&device_answers, N);

    copy_to_device(device_keys, keys, N);
    copy_to_device(device_first_child, first_child, N);
    copy_to_device(device_second_child, second_child, N);

    copy_to_device(device_search_keys, array, N);

    int threadsPerBlock = 256;
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;


    start = clock();
    cuda_flatten_tree_search<<<blocksPerGrid, threadsPerBlock>>>(N, device_search_keys, device_keys, device_first_child, device_second_child, device_answers);

    cudaMemcpy(answers_flattened_tree_cuda, device_answers, sizeof(int) * N, cudaMemcpyDeviceToHost);
    finish = clock();

    printf("searching in flattened tree data structure in CUDA took %d miliseconds\n", (finish - start) * 1000 / CLOCKS_PER_SEC);

    if (N <= 20) print_array(answers_flattened_tree_cuda, N);


    int score = 0;
    for (idx = 0; idx < N; idx++) {
        if (answers_tree[idx] == answers_flattened_tree[idx] && answers_flattened_tree[idx] == answers_flattened_tree_cuda[idx]) {
            score++;
        }
    }
    if (score == N) {
        printf("Everything is fine!\n");
    }


    cudaFree(device_keys);
    cudaFree(device_search_keys);
    cudaFree(device_first_child);
    cudaFree(device_second_child);
    cudaFree(device_answers);
    return 0;
}