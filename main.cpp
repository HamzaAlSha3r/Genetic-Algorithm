#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#define FAST cout.tie(0);cin.tie(0);ios::sync_with_stdio(0);
#define int long long
#define mm(x, v) memset(x, v, sizeof(x));
#define debug(x, i) cout << "when #i = " << i << " " << '[' << #x << " is: " << x <<  "] " << endl;
#define S second
#define F first
#define GG cout << "GGG\n";
#define MD 1000000007
#define SizeOfPopulation 500
#define MAX_NUMBER 2147483647

using namespace std;
using namespace __gnu_pbds;

template<class T> using ordered_set =tree<T, null_type, less<T>, rb_tree_tag,tree_order_statistics_node_update> ;

typedef long double ld;
typedef long long ll;


map<int,int> OriginalFreqOfTheTasks;
set<int> JobsUniquely;
vector<vector<pair<int, int>>> jobs;
vector<vector<pair<int, int>>> queueOfTasks ;
int MachinesCount;


vector<vector<pair<int, int>>> ReadingFromFile(const string& FileName) {
    ifstream file(FileName); // Open the file with the job data
    if (!file.is_open()) {
        cerr << "Failed to open file." << endl;
        exit(1);
    }
    int number_of_machines;
    file >> number_of_machines; // Read the number of machines
    MachinesCount = number_of_machines ;
    vector<vector<pair<int, int>>> jobs; // Stores all jobs
    string line; // To read each line after the first line
    // Read and discard the end of line character after reading number of machines
    getline(file, line);
    while (getline(file, line)) { // Read each job block
        getline(file, line); // Read the job details line
        istringstream iss(line);
        vector<pair<int, int>> job; // Stores tasks for current job
        string task;
        while (getline(iss, task, ',')) { // Read each task within a job
            istringstream taskStream(task);
            int machine_number, time;
            taskStream >> machine_number >> time; // Read machine number and time
            job.emplace_back(machine_number, time); // Store task in job
        }

        jobs.push_back(job); // Store job in jobs
    }
    file.close(); // Close the file
    return jobs;
}

double fitnessValue(const vector<int>& chromosome) {
    double totalTime = 0 ;
    double maxJobTime = 0 ;
    int curSubTask[jobs.size()] ;
    long long exiting[jobs.size()] ;
    for (int i = 0 ; i < jobs.size(); i++) curSubTask[i] = 0 ,exiting[i] = 0;
    int temp = MachinesCount ;
    queueOfTasks.clear() ;
    for (int i = 0 ; i < temp ; i ++ ){
        vector<pair<int,int>> temp1 ;
        queueOfTasks.push_back(temp1) ;
    }
    long long nextReadyTime[MachinesCount] ;
    for (int i = 0 ; i < MachinesCount ; i ++ ) nextReadyTime[i] = 0 ;
    for (int i = 0 ; i < chromosome.size() ; i ++ ){
        int curJob = chromosome[i] ;
        pair<int,int> task = jobs[curJob - 1][curSubTask[curJob - 1]] ;
        curSubTask[curJob - 1]++;
        int machineNumber = task.first ;
        long long emptyMachine = max(0LL , exiting[curJob - 1] - nextReadyTime[machineNumber - 1]) ;
        if (emptyMachine) queueOfTasks[machineNumber - 1].push_back({-1, emptyMachine}) ;
        queueOfTasks[machineNumber - 1].push_back({curJob, task.second}) ;
        totalTime -= (1.0/(double)jobs.size()) * exiting[curJob - 1] ;
        exiting[curJob - 1] = emptyMachine + task.second + nextReadyTime[machineNumber - 1] ;
        nextReadyTime[machineNumber - 1] = exiting[curJob - 1] ;
        totalTime += (1.0/(double)jobs.size()) * exiting[curJob - 1] ;
        if (exiting[curJob - 1] > maxJobTime) maxJobTime = exiting[curJob - 1];
    }
    return (1.0/maxJobTime) * 0.20 + (1.0/totalTime) * 0.80 ;
}


struct operatorOfChrom {
    bool operator()(const vector<int>& a, const vector<int>& b) {
        return fitnessValue(a) < fitnessValue(b) ;
    }
};

void generateInitialPopulation( priority_queue<vector<int>, vector<vector<int>>, operatorOfChrom> &population ){
    if (jobs.size() == 0 ) cout << "there is no jobs to create initial population for ." << endl ;
    vector<int> curPopulation;
    for (int i = 0 ; i < jobs.size() ; i ++ ){
        int SubTasksNum = jobs[i].size() ;
        while (SubTasksNum--) curPopulation.push_back(i + 1) ;
    }
    for(int i = 0 ; i < SizeOfPopulation ; i ++ ){
        vector <int> curChrom ;
        random_shuffle(curPopulation.begin(), curPopulation.end()) ;
        for (int j = 0; j < curPopulation.size(); j++) curChrom.push_back(curPopulation[j]) ;
        population.push(curChrom) ;
    }
}

void addRandomly(vector<vector<int>>& newPopulationToGenerate, vector<vector<int>> & curPopulation){
    vector<int> chromosome = curPopulation[0] ;
    for (int i = 0 ; i < (int)(SizeOfPopulation/10) ; i ++){
        vector<int> temp;
        random_shuffle(chromosome.begin(), chromosome.end()) ;
        for (int j = 0 ; j < chromosome.size(); j++) temp.push_back(chromosome[i]) ;
        newPopulationToGenerate.push_back(temp) ;
    }
}

int getRandomNumber(int N) { return (rand() % N) + 1;}

long double getRandomProbability(){
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dist(0.0, 1);
    return dist(gen);
}


vector<int> turnIntoValid(vector<int> chromosome){
    map<int,int> freqOfTasksInChrom;
    for(auto it : chromosome) freqOfTasksInChrom[it]++;
    vector<int> needMore, hasBonus;
    map<int,vector<int>> indicesOfSubTasks;
    for(int i = 0;i < chromosome.size(); i++) indicesOfSubTasks[chromosome[i]].push_back(i);
    for(auto it : JobsUniquely) {
        if(freqOfTasksInChrom[it] < OriginalFreqOfTheTasks[it]) needMore.push_back(it);
        else if(freqOfTasksInChrom[it] > OriginalFreqOfTheTasks[it]) hasBonus.push_back(it);
    }
    while(needMore.size() || hasBonus.size()){
        int toChooseFromNeed = getRandomNumber(needMore.size()) - 1; // index of needMore task
        int toChooseFromOverload = getRandomNumber(hasBonus.size()) - 1; // index of hasBonus task
        int tarIndexOfOverLoad = getRandomNumber(indicesOfSubTasks[hasBonus[toChooseFromOverload]].size()) - 1;
        chromosome[indicesOfSubTasks[hasBonus[toChooseFromOverload]][tarIndexOfOverLoad]] = needMore[toChooseFromNeed];
        indicesOfSubTasks[hasBonus[toChooseFromOverload]].erase(indicesOfSubTasks[hasBonus[toChooseFromOverload]].begin() + tarIndexOfOverLoad);
        freqOfTasksInChrom[needMore[toChooseFromNeed]]++;
        freqOfTasksInChrom[hasBonus[toChooseFromOverload]]--;
        if(freqOfTasksInChrom[needMore[toChooseFromNeed]] == OriginalFreqOfTheTasks[needMore[toChooseFromNeed]]) needMore.erase(needMore.begin() + toChooseFromNeed);
        if(freqOfTasksInChrom[hasBonus[toChooseFromOverload]] == OriginalFreqOfTheTasks[hasBonus[toChooseFromOverload]]) hasBonus.erase(hasBonus.begin() + toChooseFromOverload);
    }
    return chromosome;
}

vector<int> decideHalfsToCross(vector<int> &firstChrom, vector<int> &secondChrom, int toCutFrom, bool flag){
    vector<int> newChrom;
    for(int i = 0;i < firstChrom.size(); i++){
        if(i < toCutFrom){
            if(flag) newChrom.push_back(firstChrom[i]);
            else newChrom.push_back(secondChrom[i]);
        }
        else {
            if(flag) newChrom.push_back(secondChrom[i]);
            else newChrom.push_back(firstChrom[i]);
        }
    }
    return newChrom;
}

vector<int> crossOver(vector<int> &f, vector<int> &s){
    vector<int> ToReturnFunc;
    int mx = max(fitnessValue(f), fitnessValue(s));
    if(mx == fitnessValue(f)) for(auto it : f) ToReturnFunc.push_back(it);
    else for(auto it : s) ToReturnFunc.push_back(it);
    int idx = getRandomNumber(f.size());
    // when flag = 0
    vector<int> chromosome = decideHalfsToCross(f, s, idx, 0);
    chromosome = turnIntoValid(chromosome);
    if(fitnessValue(chromosome)> mx){
        ToReturnFunc.clear();
        for(auto it : chromosome) ToReturnFunc.push_back(it) ;
        mx = fitnessValue(chromosome);
    }
    // try when flag = 1
    chromosome = decideHalfsToCross(f, s, idx, 1);
    chromosome = turnIntoValid(chromosome);
    if(fitnessValue(chromosome)> mx){
        ToReturnFunc.clear();
        for(auto it : chromosome) ToReturnFunc.push_back(it) ;
        mx = fitnessValue(chromosome);
    }

    return ToReturnFunc;
}

vector<double> generateExponentialDistribution(int sz, double lambda, bool fromBack) {
    vector<double> distribution(sz);
    if(fromBack) for (int i = 0;i < sz; i++) distribution[i] = exp(-lambda * i);
    else for (int i = 0;i < sz; i++) distribution[i] = exp(lambda * i);
    double sum = accumulate(distribution.begin(), distribution.end(), 0.0);
    for(double &val : distribution) val /= sum;
    return distribution;
}

template <typename T> vector<T> selectElements(const vector<T> &previousGeneration, int N, bool fromBack) {
    vector<T> ToReturnFunc;
    vector<double> distribution = generateExponentialDistribution((int)previousGeneration.size(), 0.5, fromBack);
    bool taken[previousGeneration.size() + 1];
    for (int i = 0 ; i < previousGeneration.size() ; i++) taken[i] = false;
    default_random_engine generator(random_device{}());
    discrete_distribution<int> dist(distribution.begin(), distribution.end());
    for (int i = 0; i < N; ++i) {
        int index = dist(generator);
        while(taken[index]){
            index++;
            index %= (previousGeneration.size());
        }
        vector<int> tmpToAdd;
        for(auto it : previousGeneration[index]) tmpToAdd.push_back(it);
        ToReturnFunc.push_back(tmpToAdd);
        taken[index] = true;
    }
    return ToReturnFunc;
}

vector <vector<int>> toChooseRandomly(vector <vector<int>> &previousGeneration, double percent, bool from_back){
    vector<vector<int>> ToReturnFunc = selectElements(previousGeneration, (int)previousGeneration.size() * percent, from_back);
    return ToReturnFunc;
}

void mutate(vector<int> &Chrom){
    int times = getRandomNumber(max(1LL, (long long)Chrom.size() / 2LL - 1LL)) - 1;
    for(int i = 0;i < times; i++){
        int l = 0, r = 0;
        while(l == r || Chrom[l] == Chrom[r]) l = getRandomNumber(Chrom.size()) - 1, r = getRandomNumber(Chrom.size()) - 1;
        if(getRandomProbability() >= getRandomProbability()) swap(Chrom[l], Chrom[r]);
    }
}

void chooseToCrossOver(vector<vector<int>> &newPopulationToGenerate, vector<vector<int>> &Population, long double prob){
    vector<vector<int>> ToReturnFunc;
    for(int i = 1;i < ((int)Population.size() * 0.15); i++) ToReturnFunc.push_back(crossOver(Population[i], Population[i - 1]));
    ToReturnFunc.push_back(crossOver(Population[0], Population.back()));
    for(int i = 0;i < ((int)Population.size() * 0.85); i++){
        int l = 0, r = 0;
        while(l == r) l = getRandomNumber(Population.size()) - 1, r = getRandomNumber(Population.size()) - 1;
        ToReturnFunc.push_back(crossOver(Population[l], Population[r]));
    }
    vector<vector<int>> toAddToTheFunction = toChooseRandomly(ToReturnFunc, prob, 0);
    for(int i = 0;i < toAddToTheFunction.size(); i++) {
        vector<int> tmp;
        for(auto it : toAddToTheFunction[i]) tmp.push_back(it);
        newPopulationToGenerate.push_back(tmp);
    }
}

void chooseToMutation(vector<vector<int>> &newPopulationToGenerate, vector<vector<int>> &Population, long double prob){
    vector<vector<int>> ToReturnFunc = toChooseRandomly(Population, prob, 1);
    for(int i = 0;i < ToReturnFunc.size(); i++) {
        mutate(ToReturnFunc[i]);
        vector<int> tmp;
        for(auto it : ToReturnFunc[i]) tmp.push_back(it);
        newPopulationToGenerate.push_back(tmp);
    }
}

void chooseToReproduce(vector<vector<int>> &newPopulationToGenerate, vector<vector<int>> &Population, long double prob){
    vector<vector<int>> ToReturnFunc = toChooseRandomly(Population, prob, 0);
    for(int i = 0;i < (int)(0.1 * Population.size()); i++) {
        vector<int> tmp;
        for(auto it : Population[i]) tmp.push_back(it);
        newPopulationToGenerate.push_back(tmp);
    }
}

void addRandomly(vector<vector<int>> & newPopulationToGenerate, vector<vector<int>> & curPopulation, long double percent){
    vector<int> chromosome = curPopulation[0] ;
    for (int i = 0 ; i < (int)(SizeOfPopulation * percent) ; i++){
        vector<int> temp;
        random_shuffle(chromosome.begin(), chromosome.end()) ;
        for (int j = 0 ; j < chromosome.size(); j ++ ) temp.push_back(chromosome[j]) ;
        newPopulationToGenerate.push_back(temp);
    }
}


vector<vector<int>> fillVector(priority_queue<vector<int>, vector<vector<int>>, operatorOfChrom>& population) {
    vector<vector<int>> ToReturnFunc ;
    priority_queue<vector<int>, vector<vector<int>>, operatorOfChrom> temp ;
    while (!population.empty()) {
        vector<int> chromosome = population.top();
        population.pop();
        temp.push(chromosome);
        ToReturnFunc.push_back(chromosome);
    }
    while (!temp.empty()) {
        population.push(temp.top());
        temp.pop();
    }
    return ToReturnFunc;
}

void currentProcess(priority_queue<vector<int>, vector<vector<int>>, operatorOfChrom> &lastPopulation){
    vector<vector<int>> newPopulationToGenerate, Population = fillVector(lastPopulation);
    chooseToCrossOver(newPopulationToGenerate, Population, 0.60);
    chooseToMutation(newPopulationToGenerate, Population, 0.15);
    addRandomly(newPopulationToGenerate, Population, 0.15);
    chooseToReproduce(newPopulationToGenerate, Population, 0.25);
    while(lastPopulation.size()) lastPopulation.pop();
    for(int i = 0;i < newPopulationToGenerate.size(); i++) {
        vector<int> tmp;
        for(int j = 0;j < newPopulationToGenerate[i].size(); j++) tmp.push_back(newPopulationToGenerate[i][j]);
        lastPopulation.push(tmp);
    }
}

pair<double, double> calcFit(const vector<int>& chromosome) {
    double totalTime = 0 ;
    double maxJobTime = 0 ;
    int curSubTask[jobs.size()] ;
    long long exiting[jobs.size()] ;
    for (int i = 0 ; i < jobs.size(); i++) curSubTask[i] = 0 ,exiting[i] = 0;
    int temp = MachinesCount ;
    queueOfTasks.clear() ;
    for (int i = 0 ; i < temp ; i ++ ){
        vector<pair<int,int>> temp1 ;
        queueOfTasks.push_back(temp1) ;
    }
    long long nextReadyTime[MachinesCount] ;
    for (int i = 0 ; i < MachinesCount ; i ++ ) nextReadyTime[i] = 0 ;
    for (int i = 0 ; i < chromosome.size() ; i ++ ){
        int curJob = chromosome[i] ;
        pair<int,int> task = jobs[curJob - 1][curSubTask[curJob - 1]] ;
        curSubTask[curJob - 1]++;
        int machineNumber = task.first ;
        long long emptyMachine = max(0LL , exiting[curJob - 1] - nextReadyTime[machineNumber - 1]) ;
        if (emptyMachine) queueOfTasks[machineNumber - 1].push_back({-1, emptyMachine}) ;
        queueOfTasks[machineNumber - 1].push_back({curJob, task.second}) ;
        totalTime -= (1.0/(double)jobs.size()) * exiting[curJob - 1] ;
        exiting[curJob - 1] = emptyMachine + task.second + nextReadyTime[machineNumber - 1] ;
        nextReadyTime[machineNumber - 1] = exiting[curJob - 1] ;
        totalTime += (1.0/(double)jobs.size()) * exiting[curJob - 1] ;
        if (exiting[curJob - 1] > maxJobTime) maxJobTime = exiting[curJob - 1];
    }
    return {maxJobTime, totalTime};
}


void Print(priority_queue<vector<int>, vector<vector<int>>, operatorOfChrom> &population){
    for (int i = 0 ; i < queueOfTasks.size() ; i ++ ){
        cout << "######################## Machine: " << i + 1 << " ########################\n" ;
        int total = 0;
        cout << total << " | ";
        for (int j = 0 ; j < queueOfTasks[i].size() ; j ++ ){
            total += queueOfTasks[i][j].second;
            if(queueOfTasks[i][j].first != -1) cout << "Job" << queueOfTasks[i][j].first << " | " << total << " | " ;
            else cout << "elapsedTime" << " | " << total << " | " ;
            if(j % 10 == 0 && j) cout << "\n";
        }
        cout << "\n\n";
    }
    pair<double, double> bestFit = calcFit(population.top());
    cout << "Best Fitness Of All Time: " << fitnessValue(population.top()) << "\n";
    cout << "MaximumJobTime: " << bestFit.first << "\n";
    cout << "TotalTime: " << bestFit.second << "\n";
}

int32_t main() {
    jobs = ReadingFromFile("input.txt") ;
    priority_queue<vector<int>, vector<vector<int>>, operatorOfChrom> population ;
    for (int i = 0 ; i < jobs.size() ; i++){
        OriginalFreqOfTheTasks[i + 1] = jobs[i].size();
        JobsUniquely.insert(i + 1);
    }
    generateInitialPopulation(population);
    int Times = 10, cnt = 1;
    while(Times--) currentProcess(population);
    Print(population);
    return 0;
}
