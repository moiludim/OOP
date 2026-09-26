#include <cmath>
#include "spec_func.hpp"
#include <cstdlib>
#include <ctime>
#include <iostream>
#include "distributions.hpp"

int main(){
    
    srand(static_cast<unsigned int>(time(0)));
    
    for( int i = 0; i<100; i++){
        std::cout<<simulate_mixture(2,4,0.4)<<" ";
    }
    
    tests2();
}