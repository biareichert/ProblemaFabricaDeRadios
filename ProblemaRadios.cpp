#include <bits/stdc++.h>
#include <vector>
#include <math.h>
#include <algorithm>

/* Anotações

    Autor: Geremias Corrêa, Beatriz Reichert
    Discplina: OCEV
    Semestre: 2021/1

    Exemplo de execução: g++ -Wall -O3 -o exec nomeCodigo.cpp
        ./exec <COD> <RUN> <GEN> <POP> <DIM>

    EXECUÇÃO PADRÃO PARA ESTE PROBLEMA: ./exec 0 1 1 10 10 -- EEXECUÇÃO PADRÃO PARA TESTE
        O tamanho do cromossomo PRECISA ser 10 para o problema do exercício 4, dado a restrição 24 e 32, respectivas, para as linhas S e L
        Sendo assm, fica um cromossomo de tam 10, como 5 bits da esquerda para linha S, 5 bits da direta para linha L
        O tipo PRECISA ser 0
        O '1 1 10' é variável, setamos os 1s pois não está aplicável n execuções nem m gerações, porém o tamanhoda poupulação pode variar (no exemplo está 5)
*/

using namespace std;
using std::find;

// Realiza uma distribuição uniforme de valores entre 0 a 1 (0.00 a 1.00).
random_device randomm;
mt19937 engine{randomm()};
uniform_real_distribution<double> distribuicaoBool{0.0, 1.0};
uniform_real_distribution<double> distribuicaoInt{-5.0, 10.0};
uniform_real_distribution<double> distribuicaoReal{-10.0, 10.0};

int crossoverType = 1; //1 ´crossover de 1 corte, 2 é crossover de 2 cortes, 3 é uniform crossver (50%)

void printPopulation(vector<vector<double>> pop, int popSize, int chromosomeSize, string message){
    cout << endl << message << endl;
    for (int i = 0; i < popSize; i++){
        cout << "Cr.["<< i << "]: [";
        for (int j = 0; j < chromosomeSize; j++){
            if (j == chromosomeSize - 1) cout << pop[i][j];
            else  cout << pop[i][j] << ", ";
        }
        cout << "]" << endl;
    }
    cout << endl;
}

void printFellow(vector<double> fellow, int popSize, int chromosomeSize, string message){
    cout << endl << message << endl;
    cout << "Cr.: [";
    for (int j = 0; j < chromosomeSize; j++){
        if (j == chromosomeSize - 1) cout << fellow[j];
        else  cout << fellow[j] << ", ";
    }
    cout << "]" << endl;
}

vector<vector<double>> gerarConfiguracaoAleatoria(int codificationType, int executionsNumber, int numberGenerations, int populationSize, int chromosomeSize){
    vector<double> vetor_auxiliar;
    vector<vector<double>> vetor_final;
    uniform_real_distribution<double> distribuicaoPerm{0.0, (double) chromosomeSize};
    int cont;

    switch (codificationType){
    case 0: //binária
        for(int i = 0; i < populationSize; i++){
            for(int j = 0; j < chromosomeSize; j++){
                double aleatorio = distribuicaoBool(engine);
                if(aleatorio > 0.5)
                    vetor_auxiliar.push_back(1);
                else
                    vetor_auxiliar.push_back(0);
            }
            vetor_final.push_back(vetor_auxiliar);
            vetor_auxiliar.clear();
        }
        break;
    case 1: //int 
        for(int i = 0; i < populationSize; i++){
            for(int j = 0; j < chromosomeSize; j++){
                int aleatorio = (int) distribuicaoInt(engine);
                vetor_auxiliar.push_back(aleatorio);
            }
            vetor_final.push_back(vetor_auxiliar);
            vetor_auxiliar.clear();
        }
        break;
    case 2: //int-perm
        for(int i = 0; i < populationSize; i++){
            for(int j = 0; j < chromosomeSize; j++){
                cont = 0;
                while(1) {
                    int aleatorio = (int) distribuicaoPerm(engine);
                    for (int k : vetor_auxiliar){
                        if (aleatorio == k){
                            cont++;
                            break;
                        }
                    }
                    if (cont == 0){
                         vetor_auxiliar.push_back(aleatorio);
                         break;
                         
                    }
                    cont = 0;
                }
            }
            vetor_final.push_back(vetor_auxiliar);
            vetor_auxiliar.clear();
        }
        break;
    case 3: //real
        for(int i = 0; i < populationSize; i++){
            for(int j = 0; j < chromosomeSize; j++){
                double aleatorio = distribuicaoReal(engine);
                vetor_auxiliar.push_back(aleatorio);
            }
            vetor_final.push_back(vetor_auxiliar);
            vetor_auxiliar.clear();
        }
        break;

    default:
        return vetor_final;
    }

    return vetor_final;
}

vector<int> calc_bin(vector<double> pop, int chromosomeSize){

    int sum_s = 0;
    int sum_l = 0;
    int bits_s = 5;
    int bits_l = 5;
    int aux_s = bits_s - 1;
    int aux_l = bits_l - 1;
    vector<int> val_s_l;

    for(int i = 0; i < chromosomeSize - bits_s; i++){ //somando o número de funcionários da linha S, posição 0 a 4 (5 bits)
        if(pop[i] == 1) sum_s = sum_s + pow(2, aux_s);
        aux_s -= 1;
    }

    // Normalização linha S:
    sum_s = ceil(0 + ((24 - 0) / (pow(2, bits_s) - 1)) * sum_s);

    for(int i = chromosomeSize - bits_l; i < chromosomeSize; i++){ //somando o número de funcionários da linha L, posição 5 a 10 (6 bits)
        if(pop[i] == 1) sum_l = sum_l + pow(2, aux_l);
        aux_l -= 1;   
    }

    // Normalização linha L:
    sum_l = ceil(0 + ((16 - 0) / (pow(2, bits_l) - 1)) * sum_l);
    
    val_s_l.push_back(sum_s);
    val_s_l.push_back(sum_l);

    return val_s_l;
}


float fitness(vector<double> pop, int chromosomeSize){
        
    vector<int> val_s_l = calc_bin(pop, chromosomeSize);

    //cout << endl << "Máquinas S = " << val_s_l[0] << "; Máquinas L = " << val_s_l[1] << "; Func. utilizados: " << val_s_l[0] + 2 * val_s_l[1] << endl;
    //if (val_s_l[0] > 24) cout << "Solução inválida pois excede a restrição da quantidade limite de 24 funcionários na linha standard." << endl;
    //if (val_s_l[1] > 16) cout << "Solução inválida pois excede a restrição da quantidade limite de 32 funcionários na linha luxo." << endl;
    
    float fit = ((float(val_s_l[0]) * 30 + float(val_s_l[1]) * 40) / 1360) + (-1) * fmax(0, (float(val_s_l[0]) + 2 * float(val_s_l[1]) - 40) / 16); // fit = FOn + r * Hn
    
    return fit;
}

float fitnessVerbose(vector<double> fellow, int chromosomeSize){
        
    vector<int> val_s_l = calc_bin(fellow, chromosomeSize);

    cout << "Máquinas Standard = " << val_s_l[0] << "; Máquinas Luxe = " << val_s_l[1] << "; Total funcionários: " << val_s_l[0] + 2 * val_s_l[1] << endl;
    if (val_s_l[0] > 24) cout << "Solução inválida pois excede a restrição da quantidade limite de 24 funcionários na linha standard." << endl;
    if (val_s_l[1] > 16) cout << "Solução inválida pois excede a restrição da quantidade limite de 32 funcionários na linha luxo." << endl;
    if (val_s_l[0] + 2 * val_s_l[1] > 40) cout << "Solução inválida pois excede a restrição de 40 funcionários totais." << endl;
    
    float fit = ((float(val_s_l[0]) * 30 + float(val_s_l[1]) * 40) / 1360) + (-1) * fmax(0, (float(val_s_l[0]) + 2 * float(val_s_l[1]) - 40) / 16); // fit = FOn + r * Hn
    
    return fit;
}

vector<vector<double>> crossOverOnePoint(vector<double> fellow1, vector<double> fellow2, int chromosomeSize){

    uniform_real_distribution<double> distribuicaoInt{0.0, double(chromosomeSize)}; // não pode pegar a última posição pro 2 point, pois se não o segundo corte não tem quem pegar
    int cutIndex = distribuicaoInt(engine); //pegará até ele, incluindo ele, para lado esquerdo do vetor.

    vector<double> newFellow1;
    vector<double> newFellow2;
    vector<vector<double>> newFellows;

    for (int i = 0; i < chromosomeSize; i++){
        if (i <= cutIndex){
            newFellow1.push_back(fellow1[i]);
            newFellow2.push_back(fellow2[i]);
        }else {
            newFellow1.push_back(fellow2[i]);
            newFellow2.push_back(fellow1[i]);
        }
    }

    newFellows.push_back(newFellow1);
    newFellows.push_back(newFellow2);

    return newFellows;
}

vector<vector<double>> crossOverTwoPoint(vector<double> fellow1, vector<double> fellow2, int chromosomeSize){

    uniform_real_distribution<double> distribuicaoInt{0.0, double(chromosomeSize) - 1}; // não pode pegar a última posição, pois se não o segundo corte não tem quem pegar em alguns casos
    int cutIndex = distribuicaoInt(engine); //pegará até ele, incluindo ele, para lado esquerdo do vetor.
    uniform_real_distribution<double> distribuicaoIntIndex2{double(cutIndex) + 1, double(chromosomeSize)};
    int cutIndex2 = distribuicaoIntIndex2(engine); //pegará até ele, incluindo ele, para lado esquerdo do vetor.

    vector<double> newFellow1;
    vector<double> newFellow2;
    vector<vector<double>> newFellows;

    for (int i = 0; i < chromosomeSize; i++){
        if (i <= cutIndex){
            newFellow1.push_back(fellow1[i]);
            newFellow2.push_back(fellow2[i]);
        } else if (i <= cutIndex2){
            newFellow1.push_back(fellow2[i]);
            newFellow2.push_back(fellow1[i]);
        }
        else {
            newFellow1.push_back(fellow1[i]);
            newFellow2.push_back(fellow2[i]);
        }
    }

    newFellows.push_back(newFellow1);
    newFellows.push_back(newFellow2);

    return newFellows;
}  

vector<vector<double>> crossOverUniform(vector<double> fellow1, vector<double> fellow2, int chromosomeSize){
    //FUNFANDO
    vector<double> newFellow1;
    vector<double> newFellow2;
    vector<vector<double>> newFellows;
    
    for (int i = 0; i < chromosomeSize; i++){
        double aleatorio = distribuicaoBool(engine);
        if (aleatorio < 0.5){
            newFellow1.push_back(fellow1[i]);
            newFellow2.push_back(fellow2[i]);
        }else {
            newFellow1.push_back(fellow2[i]);
            newFellow2.push_back(fellow1[i]);
        }
    }

    newFellows.push_back(newFellow1);
    newFellows.push_back(newFellow2);

    return newFellows;
    
}  

vector<double> mutation(vector<double> fellow, int chromosomeSize){
    vector<double> newFellow;
    //cout << "Mutation:" << endl;
    
    for (int i = 0; i < chromosomeSize; i++){
        double aleatorio = distribuicaoBool(engine);
        // Probabilidade de mutação de 5%:
        if (aleatorio <= 0.05){
            //cout << i << " troca" << endl;

            if(fellow[i] == 0){
                newFellow.push_back(1);
            }else{
                newFellow.push_back(0);
            }

        }else {
            newFellow.push_back(fellow[i]);
        }
    }

    return newFellow;
}

vector<int> melhor_pior_individuo(vector<float> fitnessPerIndiv, int populationSize){
    float val_melhor_ind = fitnessPerIndiv[0];
    float val_pior_ind = fitnessPerIndiv[0];
    int posicao_melhor = 0;
    int posicao_pior = 0;
    vector<int> individuos;

    for(int i=0; i<populationSize; i++){
        if(fitnessPerIndiv[i] > val_melhor_ind){
            posicao_melhor = i;
            val_melhor_ind = fitnessPerIndiv[i];
        }

        if(fitnessPerIndiv[i] < val_pior_ind){
            posicao_pior = i;
            val_pior_ind = fitnessPerIndiv[i];
        }
    }

    individuos.push_back(posicao_melhor);
    individuos.push_back(posicao_pior);

    return individuos;
}

vector<vector<double>> selectionRoutineShortExplained(vector<vector<double>> population, int populationSize, int chromosomeSize, vector<float> fitnessPerIndiv){
    // TO DO: aqui dentro que deve ser chamado o how many queens attacked

    vector<float> probabilities; //probabilidades entre 0 e 1, flutuante.
    float totalSumFitness = 0;
    float localFitnessProportion = 0;
    float actualSumm = 0;
    int indexChoosedToFirstPair = 0;
    int indexChoosedToSecondPair = 0;

    vector<vector<double>> generatePopAfterCrossover;
    vector<vector<double>> auxFellows;

    for (int i = 0; i< populationSize; i++){
        totalSumFitness += fitnessPerIndiv[i];
    }

    //cout << endl << "Total lucros com todos cromossomos para cálculo da proporção individual do fitness: " << totalSumFitness << endl;
    //cout << "Proporção de fitness por cromossomo (fitness relativo):" << endl;
    //for (int i = 0; i< populationSize; i++){
    //    cout << "Cr[" << i << "]: " << fitnessPerIndiv[i]/totalSumFitness << endl;
    //}

    //cout << endl << "Rotina de seleção dos novos pares de indivíduos/cromossomos" << endl;

    // Neste for se itera pela população para geração dosnovos indivíduos, a partir do crossover
    for (int i = 0; i < populationSize / 2; i++){

        //Geração de cada indivíduo a partir de 1 par de indívduos
        for (int j=0; j < 2; j++){
            double aleatorio = distribuicaoBool(engine);
            actualSumm = 0;
            for (int k=0; k<populationSize; k++){ //Itera por cada individuo para avalair sua fitness para sorteio
                localFitnessProportion = fitnessPerIndiv[k] / totalSumFitness;
                actualSumm += localFitnessProportion;
                if (aleatorio < actualSumm && j == 0){ //Se primeiro elemento do par
                    indexChoosedToFirstPair = k;
                    break;
                }
                else if(aleatorio < actualSumm && j == 1 && k == indexChoosedToFirstPair ){ //Caso segundo elemento do par e repita
                    aleatorio = distribuicaoBool(engine);
                    actualSumm = 0;
                    k = 0;
                }
                else if (aleatorio < actualSumm && j == 1 && k != indexChoosedToFirstPair){ // Caso segundo elemento do par e não é repetição
                    indexChoosedToSecondPair = k;
                    break;
                }
            }
        }
        //cout << endl << "Os indivíduos/cromossomos que dão origem a este par são: os de índice " << indexChoosedToFirstPair << " e " << indexChoosedToSecondPair << endl;

        double crossOverChance = distribuicaoBool(engine);
        
        // Chance de realizar crossover: 80%
        if (crossOverChance >= 0.2){ 
            //cout << "Em it" << i << " vai realizar crossover!" << endl;

            // Aqui escolhe qual função do crossover chamar. Podemos usar variável global setada também.
            auxFellows = crossOverUniform(population[indexChoosedToFirstPair], population[indexChoosedToSecondPair], chromosomeSize);
            //auxFellows = crossOverOnePoint(population[indexChoosedToFirstPair], population[indexChoosedToSecondPair], chromosomeSize);
            //auxFellows = crossOverTwoPoint(population[indexChoosedToFirstPair], population[indexChoosedToSecondPair], chromosomeSize);
            generatePopAfterCrossover.push_back(auxFellows[0]);
            generatePopAfterCrossover.push_back(auxFellows[1]);
            auxFellows.clear();
        }else{
            generatePopAfterCrossover.push_back(population[indexChoosedToFirstPair]);
            generatePopAfterCrossover.push_back(population[indexChoosedToSecondPair]);
            //cout << "Em it" << i << " NÃO vai realizar crossover!" << endl;
        }

    }

    // IF POPULATION SIZE % 2 == 1, ENTÃO FAZ UM GERAÇÃO EXTRA SÓ COM MUTAÇÃO, SEM CROSSOVER, SOMENTE MUTAÇÃO!:
    if(populationSize % 2 == 1){
        //Giro a roleta para elemento único, sem par, enviando-a para mutação direto.
        double aleatorio = distribuicaoBool(engine);
        actualSumm = 0;
        for (int k=0; k<populationSize; k++){ //Itera por cada individuo para avalair sua fitness para sorteio
            localFitnessProportion = fitnessPerIndiv[k] / totalSumFitness;
            actualSumm += localFitnessProportion;
            if (aleatorio < actualSumm){
                indexChoosedToFirstPair = k;
                break;
            }
        }
        generatePopAfterCrossover.push_back(population[indexChoosedToFirstPair]);
    }

    //printPopulation(generatePopAfterCrossover, populationSize, chromosomeSize, "Printando população pós-crossover:");

    //Depois de passar pelo crossover vem a mutação:
    vector<vector<double>> generatePopAfterMutation;
    vector<double> mutFellows;

    for (int i = 0; i < populationSize; i++){
    
        mutFellows = mutation(generatePopAfterCrossover[i], chromosomeSize);
        generatePopAfterMutation.push_back(mutFellows);
        mutFellows.clear();
       
     }

    printPopulation(generatePopAfterMutation, populationSize, chromosomeSize, "População final desta geração:");
    
    return generatePopAfterMutation;
}  

vector<vector<double>> selectionRoutineTorneio(vector<vector<double>> population, int populationSize, int chromosomeSize, vector<float> fitnessPerIndiv){
    // TO DO: aqui dentro que deve ser chamado o how many queens attacked

    vector<float> probabilities; //probabilidades entre 0 e 1, flutuante.
    float totalSumFitness = 0;
    float localFitnessProportion = 0;
    float actualSumm = 0;
    int indexChoosedToFirstPair = 0;
    int indexChoosedToSecondPair = 0;

    vector<vector<double>> generatePopAfterCrossover;
    vector<vector<double>> auxFellows;

    // Para método do TORNEIO especificamente;
    int k = 5;
    int aux;
    double kp = 1;
    vector<int> indexes;
    vector<float> fitness;
    vector<int> m_p_ind;
    uniform_real_distribution<double> distribuicaoInt{0.0, double(populationSize)};

    for (int i = 0; i< populationSize; i++){
        totalSumFitness += fitnessPerIndiv[i];
    }

    //cout << endl << "Total lucros com todos cromossomos para cálculo da proporção individual do fitness: " << totalSumFitness << endl;
    //cout << "Proporção de fitness por cromossomo (fitness relativo):" << endl;
    //for (int i = 0; i< populationSize; i++){
    //    cout << "Cr[" << i << "]: " << fitnessPerIndiv[i]/totalSumFitness << endl;
    //}

    //cout << endl << "Rotina de seleção dos novos pares de indivíduos/cromossomos" << endl;

    // Neste for se itera pela população para geração dosnovos indivíduos, a partir do crossover
    for (int i = 0; i < populationSize / 2; i++){

        //Geração de cada indivíduo a partir de 1 par de indívduos - crossover
        for (int j=0; j < 2; j++){

            // AQUI DENTRO FICARÁ A APLICAÇÃO DO TORNEIO!

            //Passo 1: seleciona-se k indivíduos aleatoriamente. Usualmente tem-se k = 2.
            for (int x = 0; x < k; x++){
                aux = distribuicaoInt(engine);

                indexes.push_back(aux);
            }

            for (int x = 0; x < k; x++){  //Pegando a fitness dos k valores
                fitness.push_back(fitnessPerIndiv[int(indexes[x])]);
            }

            m_p_ind = melhor_pior_individuo(fitness, k);

            if(j == 0){
                if(kp >= distribuicaoBool(engine)){//Escolhe o melhor indivíduo
                    indexChoosedToFirstPair = indexes[m_p_ind[0]];
                }else{//Escolhe o pior indivíduo
                    indexChoosedToFirstPair = indexes[m_p_ind[1]];
                }
            }else{
                if(kp >= distribuicaoBool(engine)){//Escolhe o melhor indivíduo
                    indexChoosedToSecondPair = indexes[m_p_ind[0]];
                }else{//Escolhe o pior indivíduo
                    indexChoosedToSecondPair = indexes[m_p_ind[1]];
                }
            }

            indexes.clear();
            fitness.clear();
            m_p_ind.clear();
        }
        //cout << endl << "Os indivíduos/cromossomos que dão origem a este par são: os de índice " << indexChoosedToFirstPair << " e " << indexChoosedToSecondPair << endl;

        double crossOverChance = distribuicaoBool(engine);
        
        // Chance de realizar crossover: 80%
        if (crossOverChance >= 0.2){ 
            //cout << "Em it" << i << " vai realizar crossover!" << endl;

            // Aqui escolhe qual função do crossover chamar. Podemos usar variável global setada também.
            auxFellows = crossOverUniform(population[indexChoosedToFirstPair], population[indexChoosedToSecondPair], chromosomeSize);
            generatePopAfterCrossover.push_back(auxFellows[0]);
            generatePopAfterCrossover.push_back(auxFellows[1]);
            auxFellows.clear();
        }else{
            generatePopAfterCrossover.push_back(population[indexChoosedToFirstPair]);
            generatePopAfterCrossover.push_back(population[indexChoosedToSecondPair]);
            //cout << "Em it" << i << " NÃO vai realizar crossover!" << endl;
        }

    }

    // IF POPULATION SIZE % 2 == 1, ENTÃO FAZ UM GERAÇÃO EXTRA SÓ COM MUTAÇÃO, SEM CROSSOVER, SOMENTE MUTAÇÃO!:
    if(populationSize % 2 == 1){
        //Giro a roleta para elemento único, sem par, enviando-a para mutação direto.
        double aleatorio = distribuicaoBool(engine);
        actualSumm = 0;
        for (int k=0; k<populationSize; k++){ //Itera por cada individuo para avalair sua fitness para sorteio
            localFitnessProportion = fitnessPerIndiv[k] / totalSumFitness;
            actualSumm += localFitnessProportion;
            if (aleatorio < actualSumm){
                indexChoosedToFirstPair = k;
                break;
            }
        }
        generatePopAfterCrossover.push_back(population[indexChoosedToFirstPair]);
    }

    //printPopulation(generatePopAfterCrossover, populationSize, chromosomeSize, "Printando população pós-crossover:");

    //Depois de passar pelo crossover vem a mutação:
    vector<vector<double>> generatePopAfterMutation;
    vector<double> mutFellows;

    for (int i = 0; i < populationSize; i++){
    
        mutFellows = mutation(generatePopAfterCrossover[i], chromosomeSize);
        generatePopAfterMutation.push_back(mutFellows);
        mutFellows.clear();
       
     }

    printPopulation(generatePopAfterMutation, populationSize, chromosomeSize, "População final desta geração:");
    
    return generatePopAfterMutation;
}  

void objectiveFunctionOutput(double sum, vector<double> elements, double qtdElements){
    double sdPow = 0, sd = 0;
    double mean = sum / qtdElements;

    for (int i = 0; i < qtdElements; i++){
        sdPow += pow(elements[i] - mean, 2);
    }

    sd = sqrt(sdPow / qtdElements);

    cout << "Média da função objetivo com desvio padrão: " << mean << " +- " << sd << endl;
    
    cout << "## ## ## ## ## ## ## ##" << endl;
}

int main(int argc, char const *argv[]) {

    // Partindo do princípio que quero uma config de 10 variáveis
    //int quantidade_variaveis = 10;

    if (argc != 6){
		printf("5 argumentos: <COD> <RUN> <GEN> <POP> <DIM> \n");
		return 1;
	}


    int codificationType; // (COD) - tipo de codificação 0 - bin; 1 - int; 2 - int-perm; 3 - real // FIXO PARA ESTE PROBLEMA (EXERCÍCIO 4)
    int executionsNumber; // (RUN) - número de execuções
    int numberGenerations; // (GEN) - número de gerações
    int populationSize; // (POP) - tamanho da população
    int chromosomeSize; // (DIM) - tamanho do cromossomo - FIXO PARA ESTE PROBLEMA (EXERCÍCIO 4)

    sscanf(argv[1], "%d", &codificationType);
    sscanf(argv[2], "%d", &executionsNumber);
    sscanf(argv[3], "%d", &numberGenerations);
    sscanf(argv[4], "%d", &populationSize);
    sscanf(argv[5], "%d", &chromosomeSize);

    // FIXADOS PARA ESTE PROBLEMA
    //codificationType = 0;
    //chromosomeSize = 11;

    
    
    

    // SE EXECUTIONS_NUMBER == 1, TRÁS COMO ESTOU TRAZENDO, SE > 1, TRÁS AS MÉDIAS DE CADA EXECUÇÃO, MLEHOR FITNESS E MÉDIA MÉDIA FITNESS

    vector<float> bestFitnessPerGeneration;
    vector<float> meanFitnessPerGeneration;
    vector<float> meanFitnessIndivPerExecution;
    vector<float> meanMeanFitnessPopPerExecution;
    vector<float> bestFitnessPerExec;

    vector<double> bestFellowAllExec;
    float bestFitnessAllExec = 0;

    vector<double> eachObjFunctionValueAllExec;
    double sumObjFunctionValueAllExec = 0;

    // Aqui vai um for externo para cada geração possível 
    for (int numExec = 0; numExec < executionsNumber; numExec++){

        vector<vector<double>> pop_inicial = gerarConfiguracaoAleatoria(codificationType, executionsNumber, numberGenerations, populationSize, chromosomeSize);
        vector<vector<double>> pop_atual = pop_inicial;
        vector<vector<double>> pop_atual_aux;

        vector<float> fitnessPerIndiv;

        // Dados do melhor usuário:
        vector<double> bestFellow;
        float bestFitness = 0;

        float sumFitness = 0;
        float sumBestFitnessIndiv = 0;
        float sumAverageFitnessPop = 0;

        double sumObjFunctionValuePerExec = 0;
        vector<double> eachObjFunctionValuePerExec;

        // População inicial:
        cout << endl;
        cout << "########### EXECUÇÃO " << numExec << ": ########### " << endl;
        printPopulation(pop_inicial, populationSize, chromosomeSize, "Printando população inicial:");

        for(int generation = 0; generation < numberGenerations; generation ++){
            cout << endl;
            cout << "########### GERAÇÃO " << generation << ": ########### " << endl;

            printPopulation(pop_atual, populationSize, chromosomeSize, "População do início desta geração:");
        
            for (int i = 0; i < populationSize; i++){

                // FITNESS
                float val_fit = fitness(pop_atual[i], chromosomeSize);
                sumFitness += val_fit;
                fitnessPerIndiv.push_back(val_fit);
                cout << "Lucro: " << fitnessPerIndiv[i] << " (Bruto: " << fitnessPerIndiv[i] * 1360 << ")" << endl;
                
                sumObjFunctionValuePerExec += val_fit * 1360;
                sumObjFunctionValueAllExec += val_fit * 1360;
                eachObjFunctionValuePerExec.push_back(val_fit * 1360);
                eachObjFunctionValueAllExec.push_back(val_fit * 1360);
                    
            }

            vector<int> m_p_ind = melhor_pior_individuo(fitnessPerIndiv, populationSize);

            // Caso haja um novo melhor indivíduo histórico, pegamos este novo:
            if (fitnessPerIndiv[m_p_ind[0]] > bestFitness){
                bestFellow.clear();
                for (int i = 0; i < chromosomeSize; i++){  
                    bestFellow.push_back(pop_atual[m_p_ind[0]][i]);
                }
                bestFitness = fitnessPerIndiv[m_p_ind[0]];
                cout << endl;
                printFellow(bestFellow, populationSize, chromosomeSize, "Novo melhor indivíduo: ");
                cout << "Sua fitness: " << bestFitness << endl;

                // Se o melhor dessa geração for um novo melhor gerado
                bestFitnessPerGeneration.push_back(fitnessPerIndiv[m_p_ind[0]]);
                sumBestFitnessIndiv += fitnessPerIndiv[m_p_ind[0]];

                // Se o novo melhor indivíduo é melhor que o melhor de todas execs:
                if (bestFitness > bestFitnessAllExec){
                    bestFellowAllExec.clear();
                    for (int i = 0; i < chromosomeSize; i++){  
                        bestFellowAllExec.push_back(bestFellow[i]);
                    }
                    bestFitnessAllExec = bestFitness;
                }
            }
            else{
                // Caso nessa geração não haja um novo melhor, o melhor histórico substitui o pior, para garantir sua manutenção

                for (int i = 0; i < chromosomeSize; i++){ // corrigir o cromossomo do pior indivíduo
                    pop_atual[m_p_ind[1]][i] = bestFellow[i];
                }
                sumFitness = sumFitness - fitnessPerIndiv[m_p_ind[1]] + bestFitness; // corrigir sumFitness
                fitnessPerIndiv[m_p_ind[1]] = bestFitness; // corrigir fitness do pior indivíduo

                // Se o melhor dessa geração for o pior substituído pelo melhor:
                bestFitnessPerGeneration.push_back(fitnessPerIndiv[m_p_ind[1]]);
                sumBestFitnessIndiv += fitnessPerIndiv[m_p_ind[1]];
            }

            
            meanFitnessPerGeneration.push_back(sumFitness / populationSize);
            sumAverageFitnessPop += sumFitness / populationSize;
            sumFitness = 0;

            cout << endl;
            cout << "O melhor indivíduo desta geração é o Cr.["<< m_p_ind[0] <<"] com fitness = " << fitnessPerIndiv[m_p_ind[0]]<< endl;

            //Esse selection vai precisar retorna vector<vector<double>> e este será recebido no final de cada geração para usar como geração base seguinte para próximo loop.
            //pop_atual_aux = selectionRoutineShortExplained(pop_atual, populationSize, chromosomeSize, fitnessPerIndiv);
            pop_atual_aux = selectionRoutineTorneio(pop_atual, populationSize, chromosomeSize, fitnessPerIndiv);

            pop_atual.clear();
            pop_atual = pop_atual_aux;

            pop_atual_aux.clear();
            fitnessPerIndiv.clear();
            m_p_ind.clear();
        }

        meanFitnessIndivPerExecution.push_back(sumBestFitnessIndiv / numberGenerations);
        meanMeanFitnessPopPerExecution.push_back(sumAverageFitnessPop / numberGenerations);
        bestFitnessPerExec.push_back(fitness(bestFellow, chromosomeSize));

        // EXIBIÇÃO MELHOR INDIVÍDUO ENCONTRADO NESTA EXEC:
        printFellow(bestFellow, populationSize, chromosomeSize, "### MELHOR INDIVÍDUO ENCONTRADO NESTA EXECUÇÃO: ### ");
        cout << fitnessVerbose(bestFellow, chromosomeSize) << "/1 de fitness" << endl;
        cout << "Função objetivo: " << fitness(bestFellow, chromosomeSize) * 1360 << "/1360" << endl;

        // EXIBIÇÃO DA MÉDIA E DESVIO PADRÃO DA FUNÇÃO OBJETIVO:
        cout << endl << "## ESTATÍSTICAS DA FUNÇÃO OBJETIVO DESTA EXECUÇÃO: ##" << endl;
        objectiveFunctionOutput(sumObjFunctionValuePerExec, eachObjFunctionValuePerExec, populationSize * numberGenerations);
        
    }

    // EXIBIÇÃO MELHOR INDIVÍDUO ENCONTRADO TODAS EXECS:
        cout << endl;
        cout << "################### GERAL ###################" << endl;
        cout << endl;
        
        printFellow(bestFellowAllExec, populationSize, chromosomeSize, "### MELHOR INDIVÍDUO ENCONTRADO EM TODAS AS EXECUÇÕES: ### ");
        cout << fitnessVerbose(bestFellowAllExec, chromosomeSize) << "/1 de fitness " << endl;
        cout << "Função objetivo: " << fitness(bestFellowAllExec, chromosomeSize) * 1360 << "/1360" << endl;

    // EXIBIÇÃO DA MÉDIA E DESVIO PADRÃO DA FUNÇÃO OBJETIVO:
        cout << endl << "## ESTATÍSTICAS DA FUNÇÃO OBJETIVO PARA TODAS EXECUÇÕES: ##" << endl;
        objectiveFunctionOutput(sumObjFunctionValueAllExec, eachObjFunctionValueAllExec, populationSize * numberGenerations * executionsNumber);

    // RESULTADOS FINAIS E GRAVAÇÃO DE ARQUIVO:
    ofstream fitness_file;
    string output = "RadioProblem_fitnessOutputFile_" + to_string(numberGenerations) + "gens_" +to_string(executionsNumber) +"execs.txt";
    fitness_file.open(output);

    // Caso apenas uma execução, exibe e grava em arquivo o fitness do melhor indiv dde cada geração e a média do fitness da pop de cada geração:
    if (executionsNumber == 1){
        cout << endl << "LISTA DE CONVERGÊNCIA POR GERAÇÃO GERADA EM ARQUIVO DE SAÍDA!" << endl;
        cout << "Geração  | Melhor fitn.  | Média fitness" << endl;
        fitness_file << "Geracao MelhorFitn MediaFitness\n";
        for(int generation = 0; generation < numberGenerations; generation ++){
            cout << generation+1 << "\t | " << bestFitnessPerGeneration[generation] << "\t | " << meanFitnessPerGeneration[generation] << endl; 
            fitness_file << to_string(generation+1) + " " + to_string(bestFitnessPerGeneration[generation]) + " " + to_string(meanFitnessPerGeneration[generation]) + "\n";
        }
    }
    // Caso mais de uma  execução, exibe e grava em arquivo a média do fitness de cada exec e a média da média do fitness da pop de cada execução :
    else{
        cout << endl << "LISTA DE CONVERGÊNCIA POR EXECUÇÃO GERADA EM ARQUIVO DE SAÍDA!" << endl;
        cout << "Execução | Média melhor fitn.  | Média média fitness pop | Melhor fitness exec" << endl;
        fitness_file << "Execucao MediaMelhorFitn MediaMediaFitnessPop MelhorFit\n";
        for (int exec = 0; exec < executionsNumber; exec++){
            cout << exec+1 << "\t | " << meanFitnessIndivPerExecution[exec] << "\t | " << meanMeanFitnessPopPerExecution[exec] << "\t | " << bestFitnessPerExec[exec] << endl; 
            fitness_file << to_string(exec+1) + " " + to_string(meanFitnessIndivPerExecution[exec]) + " " + to_string(meanMeanFitnessPopPerExecution[exec]) + " " + to_string(bestFitnessPerExec[exec]) + "\n";
            
        }
    }
    fitness_file.close();
}

