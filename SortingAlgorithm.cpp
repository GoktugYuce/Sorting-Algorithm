#include <iostream>
#include <cstdlib>
#include <ctime>
#include <windows.h>
#include <vector>
#include <fstream>
#include <Cmath>

using std::cout;
using std::endl;

int linsearch(int* A, int victim, int len);
int binsearch(int* A, int victim, int len);
void bubbleSort(int arr[], int n);
void merge(int array[], int const left, int const mid, int const right);
void mergeSort(int array[], int const begin, int const end);
void swap(int* a, int* b);
int partition(int arr[], int low, int high);
void quickSort(int arr[], int low, int high);
void insertionSort(int arr[], int n);
void copy(int* arr1, int* arr2, int N, int initsize);
int main() {

    srand(time(NULL));

    LARGE_INTEGER t1, t2, f, diff;

    QueryPerformanceFrequency(&f);
    cout << "Reported frequence is: " << f.QuadPart << '\n'
        << "This means that a the PC ticks every " << 1000000.0 / f.QuadPart << "us"
        << "\n\n" << "Side note: The size of(LARGE_INTEGER) is "
        << sizeof(LARGE_INTEGER) << endl;

    const int N = 10, initsize = 1000, inc = 50;
    cout << "Generating the array of " << N * initsize << " bytes" << endl;
    int* A = new int[N * initsize];
    int* B = new int[N * initsize];
    A[0] = 0;
    for (int i = 1; i < N * initsize; i++) {
        B[i] = rand();
    }
    cout << "\tDONE" << endl;


    std::vector<float> bin, lin, merg, inser, bubbl, quic;

    int R; //smoothing factor, will change for each algorithm;

    for (int ArrSize = initsize; ArrSize <= N * initsize; ArrSize += inc) {
        cout << ArrSize / initsize << endl;

        copy(B, A, N, initsize);
        QueryPerformanceCounter(&t1);//Start the stopwatch

        R = 1;
        for (int j = 0; j < R; j++) {
            insertionSort(A, ArrSize);
        }

        QueryPerformanceCounter(&t2);//Stop the stopwatch
        diff.QuadPart = (t2.QuadPart - t1.QuadPart);
        //Time difference in ticks
        inser.push_back(diff.QuadPart / float(R));

        copy(B, A, N, initsize);
        QueryPerformanceCounter(&t1);
        R = 1;
        for (int j = 0; j < R; j++) {
            mergeSort(A, 0, ArrSize - 1);
        }
        QueryPerformanceCounter(&t2);
        diff.QuadPart = (t2.QuadPart - t1.QuadPart);
        //Time difference in ticks
        merg.push_back(diff.QuadPart / float(R));
        copy(B, A, N, initsize);
        QueryPerformanceCounter(&t1);
        R = 1;
        for (int j = 0; j < R; j++) {
            bubbleSort(A, ArrSize);
        }
        QueryPerformanceCounter(&t2);
        diff.QuadPart = (t2.QuadPart - t1.QuadPart);
        //Time difference in ticks
        bubbl.push_back(diff.QuadPart / float(R));
        copy(B, A, N, initsize);
        QueryPerformanceCounter(&t1);
        R = 1;
        for (int j = 0; j < R; j++) {
            quickSort(A, 0, ArrSize - 1);
        }
        QueryPerformanceCounter(&t2);
        diff.QuadPart = (t2.QuadPart - t1.QuadPart);
        //Time difference in ticks
        quic.push_back(diff.QuadPart / float(R));


    }
    delete[] A;

    std::ofstream of("r.txt", std::ofstream::out);
    of << "size\tmerg\tinser\tquic\tbubbl\n";
    for (int i = 0; i < merg.size(); i++) {
        of << initsize + i * inc
            << '\t'
            // a , as the decimal point ducttape hack
            // replace with simply lin[i], if . is fine
            << floor(merg[i]) << ',' << floor((merg[i] - floor(merg[i])) * 1000 + 0.5)
            << '\t'
            << floor(inser[i]) << ',' << floor((inser[i] - floor(inser[i])) * 1000 + 0.5)
            << '\t'
            << floor(quic[i]) << ',' << floor((quic[i] - floor(quic[i])) * 1000 + 0.5)
            << '\t'
            << floor(bubbl[i]) << ',' << floor((bubbl[i] - floor(bubbl[i])) * 1000 + 0.5)
            << '\n';
    }
    of.close();
    
    return 0;
}

//linear search algorithm finds an index at which the given element is stored
int linsearch(int* A, int victim, int len) {
    for (volatile int i = 0; i < len; i++) {
        if (A[i] == victim)
            return i;
    }
    return -1;
}

//binary search algorithm finds an index at which the given element is stored in a sorted array
int binsearch(int* A, int victim, int len) {
    volatile int mid;
    int start = 0;
    int end = len - 1;
    while (start <= end) {
        mid = (start + end) / 2;
        if (A[mid] == victim)
            return mid;
        if (A[mid] < victim)
            start = mid + 1;
        else
            end = mid - 1;
    }
    return -1;
}

void bubbleSort(int arr[], int n)
{
    int i, j;
    for (i = 0; i < n - 1; i++)

        // Last i elements are already
        // in place
        for (j = 0; j < n - i - 1; j++)
            if (arr[j] > arr[j + 1])
                std::swap(arr[j], arr[j + 1]);
}
void merge(int array[], int const left, int const mid,
    int const right)
{
    auto const subArrayOne = mid - left + 1;
    auto const subArrayTwo = right - mid;

    // Create temp arrays
    auto* leftArray = new int[subArrayOne],
        * rightArray = new int[subArrayTwo];

    // Copy data to temp arrays leftArray[] and rightArray[]
    for (auto i = 0; i < subArrayOne; i++)
        leftArray[i] = array[left + i];
    for (auto j = 0; j < subArrayTwo; j++)
        rightArray[j] = array[mid + 1 + j];

    auto indexOfSubArrayOne
        = 0, // Initial index of first sub-array
        indexOfSubArrayTwo
        = 0; // Initial index of second sub-array
    int indexOfMergedArray
        = left; // Initial index of merged array


    while (indexOfSubArrayOne < subArrayOne
        && indexOfSubArrayTwo < subArrayTwo) {
        if (leftArray[indexOfSubArrayOne]
            <= rightArray[indexOfSubArrayTwo]) {
            array[indexOfMergedArray]
                = leftArray[indexOfSubArrayOne];
            indexOfSubArrayOne++;
        }
        else {
            array[indexOfMergedArray]
                = rightArray[indexOfSubArrayTwo];
            indexOfSubArrayTwo++;
        }
        indexOfMergedArray++;
    }

    while (indexOfSubArrayOne < subArrayOne) {
        array[indexOfMergedArray]
            = leftArray[indexOfSubArrayOne];
        indexOfSubArrayOne++;
        indexOfMergedArray++;
    }

    while (indexOfSubArrayTwo < subArrayTwo) {
        array[indexOfMergedArray]
            = rightArray[indexOfSubArrayTwo];
        indexOfSubArrayTwo++;
        indexOfMergedArray++;
    }
    delete[] leftArray;
    delete[] rightArray;
}


void mergeSort(int array[], int const begin, int const end)
{
    if (begin >= end)
        return;

    auto mid = begin + (end - begin) / 2;
    mergeSort(array, begin, mid);
    mergeSort(array, mid + 1, end);
    merge(array, begin, mid, end);
}
void swap(int* a, int* b)
{
    int t = *a;
    *a = *b;
    *b = t;
}
int partition(int arr[], int low, int high)
{
    int pivot = arr[high];
    int i
        = (low
            - 1);

    for (int j = low; j <= high - 1; j++) {

        if (arr[j] < pivot) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}
void quickSort(int arr[], int low, int high)
{
    if (low < high) {

        int pi = partition(arr, low, high);


        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}
void insertionSort(int arr[], int n)
{
    int i, key, j;
    for (i = 1; i < n; i++)
    {
        key = arr[i];
        j = i - 1;


        while (j >= 0 && arr[j] > key)
        {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}

void copy(int* arr1, int* arr2, int N, int initsize) {
    for (int i = 0;i < N * initsize;i++) {
        arr2[i] = arr1[i];
    }
}
