#ifndef COMPLICATEDCONVERTOR_H
#define COMPLICATEDCONVERTOR_H

#include <iostream>

class ComplicatedConvertor
{
public:
    ComplicatedConvertor() {}
    template<typename T, typename U>
    static void ConvertVector(std::function<U()>* f, std::vector<*T> & input, std::vector<U> & output)
    {
        output.resize(input.size());
        for (int i = 0; i < input.size(); i++)
        {
        output[i] = input[i] -> (*f);
        }
    }

protected:

private:
};

#endif // COMPLICATEDCONVERTOR_H
