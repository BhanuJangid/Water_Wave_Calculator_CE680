#include<iostream>
#include<cmath>
using namespace std;

int main(){

float g = 9.8;

float input_l;
cout<<"input the assumed L"<<endl;
cin>>input_l;

float t;
cout<<"input time :"<<endl;
cin>>t;

float d;
cout<<"input d:"<<endl;
cin>>d;

float cal_l=0;
int iteration_count = 0;

while(true){

    cal_l = ((g*pow(t,2))/(2*3.14))*(tanh((2*3.14*d/input_l)));
    float error = abs(cal_l- input_l);
    cout<<cal_l<<endl;
    iteration_count++;

    if(error<0.001) break;
    input_l = cal_l;

}

cout<<"Answer: "<<cal_l<<endl;
cout<<"no of iterations:    "<<iteration_count<<endl;


return 0;
}