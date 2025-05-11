#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>

//CREDITS: "https://github.com/portfoliocourses/c-example-code/blob/main/insertion_sort.c"
// "https://www.programiz.com/dsa/heap-sort"
// "https://www.geeksforgeeks.org/quick-sort-in-c/"

#define INSERTION "insertion"
#define HEAP "heap"
#define QUICK "quick"
#define QUANTIDADETAMANHOS 10
#define ULTIMOINDICE QUANTIDADETAMANHOS-1
#define MAXARRAY 256
#define VALOR_MAXIMO 20

int tamanhoArrays1[] = {10000,20000,30000,40000,50000,60000,70000,80000,90000,100000};
int tamanhoArrays[] = {100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000};

float arrayTemposA[] = {0,0,0,0,0,0,0,0,0,0};
float arrayTemposB[] = {0,0,0,0,0,0,0,0,0,0};
float arrayTemposC[] = {0,0,0,0,0,0,0,0,0,0};

void gerarGrafico(char *nomeAlgoritmo, char *nomeFicheiro);
void quicksort(int array[], int length);
void quicksort_recursion(int array[], int low, int high);
int partition(int array[], int low, int high);

void salvarDados(int n, char *nomeAlgoritmo) {
    char nomeFicheiro[MAXARRAY];
    snprintf(nomeFicheiro, sizeof(nomeFicheiro), "%s.txt", nomeAlgoritmo);

    FILE *arquivoTempos = fopen(nomeFicheiro, "w");
    if (!arquivoTempos) {
        printf("Erro ao abrir o arquivo!\n");
        return;
    }
    
    for (int i = 0; i < n; i++) {
        fprintf(arquivoTempos, "%d %.6f %.6f %.6f\n", tamanhoArrays[i], arrayTemposA[i], arrayTemposB[i], arrayTemposC[i]);
    }

    fclose(arquivoTempos);

    gerarGrafico(nomeAlgoritmo, nomeFicheiro);
}

void gerarGrafico(char *nomeAlgoritmo, char *nomeFicheiro) {
    FILE *gnuplot = popen("gnuplot -persistent", "w");
    if (!gnuplot) {
        printf("Erro ao abrir o GNUplot!\n");
        return;
    }

    int maximoX = tamanhoArrays[ULTIMOINDICE] + 10000;
    

    fprintf(gnuplot, "set title 'Tempo de Execucao do %sSort'\n", nomeAlgoritmo);
    fprintf(gnuplot, "set xlabel 'Tamanhos'\n");
    fprintf(gnuplot, "set ylabel 'Tempo (s)'\n");
    fprintf(gnuplot, "set key left top\n");

    float maximoA = arrayTemposA[ULTIMOINDICE];
    float maximoB = arrayTemposB[ULTIMOINDICE];
    float maximoC = arrayTemposC[ULTIMOINDICE];
    float maximoY;

    if(fmaxf(maximoA, maximoB) == maximoA && fmaxf(maximoA, maximoC)){
        maximoY = maximoA;
    } else if(fmaxf(maximoB, maximoC) == maximoB){
        maximoY = maximoB;
    } else {
        maximoY = maximoC;
    }

    fprintf(gnuplot, "set yrange [0:%.6f]\n", maximoY);
    fprintf(gnuplot, "set xrange [%d:%d]\n", tamanhoArrays[0], maximoX);
    fprintf(gnuplot, "set xtics 0,100000,%d\n", tamanhoArrays[ULTIMOINDICE]);

    if(strcmp(nomeAlgoritmo, INSERTION) == 0){
        fprintf(gnuplot, "regressaoA(x) = a1*x**2\n");
        fprintf(gnuplot, "regressaoB(x) = a2*x**2\n");
        fprintf(gnuplot, "regressaoC(x) = a3*x**2\n");

        fprintf(gnuplot, "a1 = 1e-10\n");
        fprintf(gnuplot, "a2 = 1e-10\n");
        fprintf(gnuplot, "a3 = 1e-10\n");

        fprintf(gnuplot, "fit regressaoA(x) '%s' using 1:2 via a1\n", nomeFicheiro);
        fprintf(gnuplot, "fit regressaoB(x) '%s' using 1:3 via a2\n", nomeFicheiro);
        fprintf(gnuplot, "fit regressaoC(x) '%s' using 1:4 via a3\n", nomeFicheiro);

        fprintf(gnuplot, "plot '%s' using 1:2 with points lc rgb 'red' title 'arrayA', regressaoA(x) with lines lc rgb 'red' title 'Regressao n**2 A', \\\n", nomeFicheiro);
        fprintf(gnuplot, "     '%s' using 1:3 with points lc rgb 'blue' title 'arrayB', regressaoB(x) with lines lc rgb 'blue' title 'Regressao n**2 B', \\\n", nomeFicheiro);
        fprintf(gnuplot, "     '%s' using 1:4 with points lc rgb 'green' title 'arrayC', regressaoC(x) with lines lc rgb 'green' title 'Regressao n**2 C'\n", nomeFicheiro);
    } else if(strcmp(nomeAlgoritmo, HEAP) == 0 || strcmp(nomeAlgoritmo, QUICK) == 0){
        fprintf(gnuplot, "regressaoA(x) = a1*x*log10(x)\n");
        fprintf(gnuplot, "regressaoB(x) = a2*x*log10(x)\n");
        fprintf(gnuplot, "regressaoC(x) = a3*x*log10(x)\n");

        fprintf(gnuplot, "a1 = 1e-7\n");
        fprintf(gnuplot, "a2 = 1e-7\n");
        fprintf(gnuplot, "a3 = 1e-7\n");

        fprintf(gnuplot, "fit regressaoA(x) '%s' using 1:2 via a1\n", nomeFicheiro);
        fprintf(gnuplot, "fit regressaoB(x) '%s' using 1:3 via a2\n", nomeFicheiro);
        fprintf(gnuplot, "fit regressaoC(x) '%s' using 1:4 via a3\n", nomeFicheiro);

        fprintf(gnuplot, "plot '%s' using 1:2 with points lc rgb 'red' title 'arrayA', regressaoA(x) with lines lc rgb 'red' title 'Regressao n*Log(n) A', \\\n", nomeFicheiro);
        fprintf(gnuplot, "     '%s' using 1:3 with points lc rgb 'blue' title 'arrayB', regressaoB(x) with lines lc rgb 'blue' title 'Regressao n*Log(n) B', \\\n", nomeFicheiro);
        fprintf(gnuplot, "     '%s' using 1:4 with points lc rgb 'green' title 'arrayC', regressaoC(x) with lines lc rgb 'green' title 'Regressao n*Log(n) C'\n", nomeFicheiro);
    }

    pclose(gnuplot);
}

void shuffle(int *array, int size) { //função para misturar os valores dos arrays
    for (int i = size - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }
}

//-----------------------------------------------------------------------------------
void insertionSort(int array[], int tamanhoArray){
    // percorre cada elemento do array a partir do segundo elemento
    // (o primeiro elemento é considerado já ordenado)
    for (int i = 1; i < tamanhoArray; i++)
    {
        // pega no elemento selecionado e vai comparar com os outros elementos que estão à sua esquerda e vai trocando-os caso os à esquerda sejam maiores
        int key = array[i];
        int j = i - 1;
        while (j >= 0 && array[j] > key)
        {
        array[j + 1] = array[j];
        j = j - 1;
        }
        array[j + 1] = key;
    }

}
//------------------------------------------------------------------------------------

void swap(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

void heapify(int arr[], int n, int i) { //função para meter array como max heap tree
    // encontra o maior valor entre raiz e os seus filhos
    int largest = i; //indice do pai
    int left = 2 * i + 1; //indice do filho esquerdo
    int right = 2 * i + 2; //indice do filho direito

    if (left < n && arr[left] > arr[largest]) //se o filho esquerdo for maior que o pai
        largest = left;

    if (right < n && arr[right] > arr[largest]) //se o filho direito for maior que o pai
        largest = right;

    // se a raiz não for o maior, vai trocar o maior nó com a raiz e voltar a ver se já está max heap
    if (largest != i) {
        swap(&arr[i], &arr[largest]);
        heapify(arr, n, largest);
    }
}

void heapSort(int arr[], int tamanho) {
    // constroi a árvore max heap
    for (int i = tamanho / 2 - 1; i >= 0; i--)
        heapify(arr, tamanho, i);

    // Heap sort
    for (int i = tamanho - 1; i >= 0; i--) {
        swap(&arr[0], &arr[i]); //mete a raiz no último nó da árvore para tirar 

        //faz outra vez heap tree
        heapify(arr, i, 0);
    }
}
//----------------------------------------------------------------------------------

int partition(int arr[], int low, int high) {

    int meio = (low + high) / 2;

    // Ordena os três elementos: arr[low], arr[meio], arr[high]
    if (arr[low] > arr[meio]) {
        swap(&arr[low], &arr[meio]);
    }
    if (arr[low] > arr[high]) {
        swap(&arr[low], &arr[high]);
    }
    if (arr[meio] > arr[high]) {
        swap(&arr[meio], &arr[high]);
    }

    int pivot_value = arr[meio];
    //colocar o pivot para o penúltimo faz a execução um pouco mais lenta
    swap(&arr[meio], &arr[high-1]); // Move o pivô para o final
    int i = low; 
    
    // não inclui o pivot
    for (int j = low; j < high-2; j++){
        // coloca à frente do pivot se o valor de j for menor ou igual ao pivot
        if (arr[j] <= pivot_value){
            swap(&arr[i], &arr[j]);
            i++;
        }
    }
    
    // mete o pivot no lugar certo 
    swap(&arr[i], &arr[high-1]);
    
    // valor do pivot
    return i;
}

void quickSort(int arr[], int low, int high) {
    while (low < high) {
        // insertion sort para arrays pequenos
        if (high - low + 1 < 1000) {
            insertionSort(arr + low, high - low + 1);
            return;
        }

        // Particiona o array e obtém o índice do pivô
        int pi = partition(arr, low, high);

        // Recursão para os subarrays
        if (pi - low < high - pi) {
            quickSort(arr, low, pi - 1);
            low = pi + 1; // Atualiza o limite inferior para evitar recursão
        } else {
            quickSort(arr, pi + 1, high);
            high = pi - 1; // Atualiza o limite superior para evitar recursão
        }
    }
}

//---------------------------------------------------------------------------------------------

int main(){
    int tamanho, escolha;
    printf("Introduza o tipo de algoritmo de ordenamento que pretende efetuar para os arrays:\n1 - Insertion Sort\n2 - Heap Sort\n3 - Quick Sort\nOpção -> ");
    scanf("%d", &escolha);
    printf("\n");

    int *arrayA;
    int *arrayB;
    int *arrayC;

    clock_t start, end;

    srand(time(NULL)); 

    for(int i = 0; i < QUANTIDADETAMANHOS; i++){
        if(escolha == 1){
            tamanho = tamanhoArrays1[i];
        } else {
            tamanho = tamanhoArrays[i]; //obtém o tamanho do array
        }

        arrayA = malloc(tamanho * sizeof(int));
        arrayB = malloc(tamanho * sizeof(int));
        arrayC = malloc(tamanho * sizeof(int));
        arrayA[0] = 0;

        for(int l = 1; l < tamanho; l++){
            if(rand() % 100 + 1 <= 5){
                arrayA[l] = arrayA[l - 1];
            } else {
                arrayA[l] = arrayA[l - 1] + 1; //adiciona valores ao array
            }
        }

        //--------------------------------------------------------------

        arrayB[0] = tamanho;

        for(int l = 1; l < tamanho; l++){
            if(rand() % 100 + 1 <= 5){
                arrayB[l] = arrayB[l - 1];
            } else {
                arrayB[l] = arrayB[l - 1] - 1;
            }
        }

        //------------------------------------------------------------

        arrayC[0] = 0;

        for(int l = 1; l < tamanho; l++){
            if(rand() % 100 + 1 <= 5){
                arrayC[l] = arrayC[l - 1];
            } else {
                arrayC[l] = arrayC[l - 1] + 1;
            }
        }

        shuffle(arrayC, tamanho);   

        switch (escolha){
        case 1:
            start = clock();
            insertionSort(arrayA, tamanho);
            end = clock();
            arrayTemposA[i] = (float)(end-start) / CLOCKS_PER_SEC;

            start = clock();
            insertionSort(arrayB, tamanho);
            end = clock();
            arrayTemposB[i] = (float)(end-start) / CLOCKS_PER_SEC;

            start = clock();
            insertionSort(arrayC, tamanho);
            end = clock();
            arrayTemposC[i] = (float)(end-start) / CLOCKS_PER_SEC;
            break;

        case 2:
            start = clock();
            heapSort(arrayA, tamanho);
            end = clock();
            arrayTemposA[i] = (float)(end-start) / CLOCKS_PER_SEC;

            start = clock();
            heapSort(arrayB, tamanho);
            end = clock();
            arrayTemposB[i] = (float)(end-start) / CLOCKS_PER_SEC;

            start = clock();
            heapSort(arrayC, tamanho);
            end = clock();
            arrayTemposC[i] = (float)(end-start) / CLOCKS_PER_SEC;
            break;

        case 3:
            start = clock();
            quickSort(arrayA, 0, tamanho-1);
            end = clock();
            arrayTemposA[i] = (float)(end-start) / CLOCKS_PER_SEC;

            start = clock();
            quickSort(arrayB,0, tamanho-1);
            end = clock();
            arrayTemposB[i] = (float)(end-start) / CLOCKS_PER_SEC;

            start = clock();
            quickSort(arrayC, 0, tamanho-1);
            end = clock();
            arrayTemposC[i] = (float)(end-start) / CLOCKS_PER_SEC;
            break;
        
        default:
            break;
        }

        free(arrayA);
        free(arrayB);
        free(arrayC);
    }
    switch (escolha){
    case 1:
        salvarDados(QUANTIDADETAMANHOS, INSERTION);
        break;

    case 2:
        salvarDados(QUANTIDADETAMANHOS, HEAP);
        break;
    
    case 3:
        salvarDados(QUANTIDADETAMANHOS, QUICK);
        break;

    default:
        break;
    }
    
    printf("\n");

    return 0;
}
