#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>

#define ORDER 1
#define REV_ORDER 2
#define BITS_IN_BYTE 8
#define MIN_COUNT 10
#define MAX_COUNT 10000
#define ARR_COUNT 4

int n_move, n_cmp;

void swap(double arr[], int i, int j)
{
    double tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
}

double random_num(void)
{
    double num;
    int *p = (int*)&num;
    *p = 0;
    for(int i = sizeof(double) / sizeof(int); i > 0; i--)
    {
      for(int j = sizeof(int) * BITS_IN_BYTE; j > 0; j--)
      {
         *p <<= 1;
         *p |= rand() % 2;
      }
      p++;
    }
    return num;
}

double *gen_arr(int type, int count)
{
    double *arr = malloc(sizeof(double) * count);
    if (type == ORDER)
    {
        for(int i = 0; i < count; i++)
            arr[i] = i;
        return arr;
    }
    else if (type == REV_ORDER)
    {
        for(int i = 0; i < count; i++)
            arr[i] = count - i;
        return arr;
    }
    else
    {
        for(int i = 0; i < count; i++)
            arr[i] = random_num();
        return arr;
    }
}

double *arr_cpy(double *arr, int count)
{
    double *new_arr = malloc(sizeof(double) * count);
    for(int i = 0; i < count; i++)
        new_arr[i] = arr[i];
    return new_arr;
}

int cmp(double x, double y)
{
    if (fabs(x) > fabs(y))
        return 1;
    return 0;
}

void Shell_Sort(double arr[], int count)
{
	for (int step = count / 2; step > 0; step /= 2)
        for (int i = step; i < count; i++)
			for (int j = i - step; j >= 0 && (n_cmp++, cmp(arr[j], arr[j + step])); j -= step)
            {
                swap(arr, j, j + step);
                n_move++;
            }
}

void Simple_Choose(double arr[], int size)
{
    for (int i = 0; i < size - 1; i++)
    {
        int min_i = i;
        for(int j = i + 1; j < size; j++)
        {
            if (cmp(arr[min_i], arr[j]))
                min_i = j;
            n_cmp++;
        }
        if (min_i != i)
        {
            swap(arr, i, min_i);
            n_move++;
        }
    }
}

void test(double* arr, int count)
{
    for(int i = 0; i < count - 1; i++)
        if (fabs(arr[i]) > fabs(arr[i + 1]))
        {
            printf("ERROR: %e %e\n", arr[i], arr[i + 1]);
            break;
        }
}

void print_arr(double arr[], int count)
{
    for(int i = 0; i < count; i++)
        printf("%e ", arr[i]);
    printf("\n");
}

int main(void)
{
    srand(time(NULL));
    for(int i = MIN_COUNT; i <= MAX_COUNT; i *= 10)
    {
        for(int j = 1; j <= ARR_COUNT; j++)
        {
            double *array1 = gen_arr(j, i);
            double *array2 = arr_cpy(array1, i);
            if (i == MIN_COUNT)
                print_arr(array1, MIN_COUNT);
            Shell_Sort(array1, i);
            printf("Shell: cmp = %d, move = %d, order = %d, count = %d\n", n_cmp, n_move, j, i);
            if (i == MIN_COUNT)
                print_arr(array1, MIN_COUNT);
            n_cmp = n_move = 0;
            Simple_Choose(array2, i);
            printf("Simple: cmp = %d, move = %d, order = %d, count = %d\n", n_cmp, n_move, j, i);
            if (i == MIN_COUNT)
                print_arr(array1, MIN_COUNT);
            n_cmp = n_move = 0;
            free(array1);
            free(array2);
        }
        printf("\n");
    }
    return 0;
}
