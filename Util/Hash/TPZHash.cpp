
#include "TPZHash.h"

uint32_t Hash(std::string str) {

    std::string::size_type size = str.size();

    uint8_t *arr = new uint8_t[size];

    for (unsigned int i = 0; i < size; ++i) {
		arr[i] = str.at(i);
    }

    uint32_t out;

    MurmurHash3_x86_32(arr, size, 0xB0F57EE3, &out);

	delete[] arr;

    return out;
}

template <>
int ClassIdOrHash<TPZFlopCounter>(){
    return Hash("TPZFlopCounter");
}

template <>
int ClassIdOrHash<int>(){
    return Hash("int");
}

template <>
int ClassIdOrHash<long int>(){
    return Hash("long int");
}

template <>
int ClassIdOrHash<float>(){
    return Hash("float");
}

template <>
int ClassIdOrHash<double>(){
    return Hash("double");
}

template <>
int ClassIdOrHash<long double>(){
    return Hash("long double");
}

template <>
int ClassIdOrHash<std::complex<float>>(){
    return Hash("std::complex<float>");
}

template <>
int ClassIdOrHash<std::complex<double>>(){
    return Hash("std::complex<double>");
}

template <>
int ClassIdOrHash<std::complex<long double>>(){
    return Hash("std::complex<long double>");
}