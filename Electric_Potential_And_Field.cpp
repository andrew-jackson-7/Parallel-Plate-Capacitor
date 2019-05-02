// This program numerically solves for the electrostatic potential and field
// from a parallel plate capacitor.
//
// The program reads in 3 parameters from a text file:
//
// "L" is the positive integer number of "Nplate" unit steps taken
// to create the parallel plate length.
//
// "Lbox" is the positive integer number of "Nplate" unit steps taken
// to create the square box encapsulating the parallel plate.
//
// "Nplate" is the unit length used to create the lengths of the parallel
// plates and the box encapsulating the parallel plates.
//
// A 2 dimensional lattice is created and indexed, where each point in the lattice
// is an "Nplate" distance apart. The entire lattice can be thought of as a matrix
// where each lattice point can be mapped to an indexed matrix element.
//
// Initial electrostatic potential values are assigned to all lattice points,
// 0 for empty space
// -1 or 1 for either capacitor plate.
//
// Subsequent electrostatic potential values are calculated by computing
// the average of the 4 nearest neighbor points.
//
// Boundary conditions were not specified at the encapsulating box edges,
// so I will not be calculating subsequent potential values there.
//
// In spite of this, the final approximation of the potential is still good.
//
// The squared diagonal elements are summed with each iteration. The
// previous and current iteration sums are compared and when their difference
// falls below an arbitrary set tolerance, the program will stop and write out
// the potential to a text file.
//
// Subsequently, the electrostatic field is calculated from the final returned
// potential using a finite difference formula.
//
// The electrostatic potential and field are simulated in Mathematica.

#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>

// The definition below maps each point on the 2 dimensional lattice to
// an element in a 1 dimensional array.

#define V(i,j) (Potential[((i) + Nbox) * (2 * Nbox + 1) + ((j) + Nbox)])

// Several using-declarations for the purpose of writing out and reading
// in from a text file.

// using std::cin;
using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::stringstream;
using std::ofstream;

// Function declarations

double Calc_Trace(double* Trace_Array_Vals, int Nbox);
double PercentDifference(double old_trace, double new_trace);
void Electric_Field_x(int Nbox, int Nplate, double* Potential);
void Electric_Field_y(int Nbox, int Nplate, double* Potential);

int main(int argc, char** argv)
{

    int count = 0;
    
// "dataIn" is an object of the ifstream class (Input stream class to operate on
// files).
//
// "store" is a string object that will hold a line that has been read in from
// the textfile.
//
// "dataIn.open("Params.txt")" calls a member function to open the file
// "Param.txt".
//
// The first while loop below cycles line by line through "params.txt"
// and counts the number of lines. In this case, "Params.txt" only has 3
// parameters on each line and therefore the count will be 3.
    
    string store;
    ifstream dataIn;
    dataIn.open("Params.txt");
    
    while (getline(dataIn, store))
    {
        count++;
    }
    
// "Parameters" is a dynamically declared array of double type and size "count".
//
// "dataIn.clear()" removes a flag when end of file is reached so that I can
// perform "getline(...)" again.
//
// "dataIn.seekg(0, dataIn.beg)" will seek to 0 characters from the beginning of
// the file.
//
// The while loop serves the same purpose of cycling through the file line by
// line and placing them into the Parameters array.
    
    string* Parameters = new string[count];
    dataIn.clear();
    dataIn.seekg(0, dataIn.beg);
    
    int newcount = 0;
    
    while (getline(dataIn, store))
    {
        Parameters[newcount] = store;
        newcount++;
    }
    
// The while loop below cycles through the elements of Parameters
// An object of stringstream is constructed and its contents are read into
// the numeric variable "Result". Each element of ParameterVals is then assigned
// to each "Result" from a while loop iteration.
    
    int Result;
    int ParameterVals[count];
    
    for (int a = 0; a < count; a++)
    {
        stringstream convert(Parameters[a]);
        convert >> Result;
        ParameterVals[a] = Result;
    }
    
// "L", "Lbox", and "Nplate" are the parameters read in from "Params.txt".
    
    int L = ParameterVals[0];
    int Lbox = ParameterVals[1];
    int Nplate = ParameterVals[2];
    
// "Nlplate" is half the length of a parallel plate centered at the origin.
// "Nbox" is half the length of a side of the square encapsulating the
// parallel plate.
    
    int Nlplate = L * Nplate;
    int Nbox = Lbox * Nplate;
    
// "Potential" is a dynamically declared array of double type with a size
// equal to the number of lattice points.
//
// "old_trace" is a numeric variable to hold the sum of the squared diagonal
// lattice point potential values.
//
// The nested for loops cycle through all lattice points and assign an initial
// electrostatic potential to each.
//
// I chose to assign "0" to all lattice points in space and "-1" or "1" to
// all plate points.
//
// Depending on the intial point values, the potential values will converge
// faster.
    
    double* Potential = new double[(2 * Nbox + 1) * (2 * Nbox + 1)];
    double old_trace;
    
    for (int i = -Nbox; i <= Nbox; i++)
    {
    
        for (int j = -Nbox; j <= Nbox; j++)
        {
        
            if ((i >= -Nlplate && i <= Nlplate) && (j == -Nlplate / L))
            {
            
                V(i, j) = -1;

            }
            
            else if ((i >= -Nlplate && i <= Nlplate) && (j == Nlplate / L))
            {
            
                V(i, j) = 1;
                
            }
            
            else
            {
            
                V(i, j) = 0;
                
            }
            
        }
        
    }
    
    old_trace = 0;
    
// The for loop below sums the squared diagonal values of the lattice (matrix)
// where they are assigned to "old_trace".
    
    for (int f = -Nbox + 1; f <= Nbox - 1; f++)
    {
    
        old_trace = pow(V(f, f), 2) + old_trace;
    
    }
    
// "new_trace" is a numeric variable to hold the most current squared diagonal
// elements sum.
//
// "Trace_Array_Vals" is a dynamically declared array of double type to hold
// the squared elements. The array does not include 2 corner points because
// they do not have boundary conditions assigned.
//
// "new_trace" is initialized to "3" (arbitrary value, could have been anything).
    
    double new_trace;
    double Trace_Array_Vals[2 * Nbox - 1];
    
    new_trace = 3;
    
// The while loop below compares "old_trace" and "new_trace" using a function
// called "PercentDifference" such that it keeps iterating until the returned
// difference falls below some tolerance, in this case ".001".
//
// The nested for loops recalculate the potential by taking the average of the
// 4 nearest neighbor points for all lattice points not on the edge
// of the encapsulating square (except parallel plate points, they are fixed).
//
// Also, the squared diagonal elements are stored in Trace_Array_Vals and then
// used by a function called "Calc_Trace" which returns the sum. This new sum
// is then assigned to "new_trace".
//
// dataOut5 is an ofstream object to write out the "PercentDifference" per iteration
// to "PercentDifference.txt".

    int counter = 1;
    
    ofstream dataOut5("PercentDifference.txt");
    
    while (PercentDifference(old_trace, new_trace) >= .001)
    {
        
        old_trace = new_trace;
        
        for (int i = -Nbox + 1; i <= Nbox - 1; i++)
        {
        
            Trace_Array_Vals[i + Nbox] = V(i, i) * V(i, i);
            
            for (int j = -Nbox + 1; j <= Nbox - 1; j++)
            {
            
                if ((i >= -Nlplate && i <= Nlplate) && (j == -Nlplate / L))
                {
                
                    V(i, j) = -1;
                    
                }
                
                else if ((i >= -Nlplate && i <= Nlplate) && (j == Nlplate / L))
                {
                
                    V(i, j) = 1;
                    
                }
                
                else
                {
                
                    V(i, j) = (V(i + 1, j) + V(i - 1, j) + V(i, j + 1) + V(i, j - 1)) / 4;
                    
                }
                
            }
            
        }
       
        new_trace = Calc_Trace(Trace_Array_Vals, Nbox);
        
        dataOut5 << counter << " " << old_trace << " " << new_trace << " " << PercentDifference(old_trace, new_trace) << endl;
        
        counter = counter + 1;
    
    }
    
// dataOut is an object that holds the "Potential_Values.txt" text file which
// is created on the fly.
//
// The potential values from the final iteration are written out to the file.
    
    ofstream dataOut("Potential_Values.txt");
    
    for (int i = -Nbox; i <= Nbox; i++)
    {
    
        for (int j = -Nbox; j <= Nbox; j++)
        {
        
            dataOut << i << " " << j << " " << V(i, j) << endl;
            
        }
        
    }
    
// "dataOut.close()" is a member function used to close the stream to the
// output file so it can't be written to without being opened again.
    
    dataOut.close();
    
// "Electric_Field_x" and "Electric_Field_y" are functions that subsequently
// calculate and write out the x and y components of the electrostatic field
// to text files (will create a text file on the fly).
    
    Electric_Field_x(Nbox, Nplate, Potential);
    Electric_Field_y(Nbox, Nplate, Potential);
    
    ofstream dataOut4("VectorPlot.txt");
    
// The x and y components of the electrostatic field are computed for each lattice
// site using the finite difference formula. Because the length between every site
// is 1, the formula is just the difference in potential between neighboring sites.
// With the components, the field is simulated in Mathematica.
    
    for (int i = -Nbox; i <= Nbox; i++)
    {
    
        for (int j = -Nbox; j <= Nbox; j++)
        {
        
            dataOut4 << i << " " << j << " " << -(V(i + 1, j) - V(i, j)) << " " << -(V(i, j + 1) - V(i, j)) << endl;
        }
        
    }
    
    return EXIT_SUCCESS;
    
}

// "Calc_Trace" computes and returns the sum of the squared diagonal elements from the lattice.
//
// Parameters:
//
// "Trace_Array_Vals" is a pointer to an array holding the non-edge squared diagonal
// elements.
//
// "Nbox" is an integer value corresponding to half the total length of one side of the square
// encapsulating the parallel plates.
//
// Return Value:
//
// "trace_count" is a numeric value holding the sum of squared diagonal elements and is returned
// by the function.

double Calc_Trace(double* Trace_Array_Vals, int Nbox)
{

    double trace_count = 0.0;
    
    for (int b = 0; b < 2 * Nbox - 2; b++)
    {
       
        trace_count = trace_count + Trace_Array_Vals[b];
        
    }
    
    return trace_count;
    
}

// "PercentDifference" calculates the decimal difference as a gauge to see how much the
// previous sum of the squared diagonal elements ("old_trace") compares to the new sum ("new_trace").
//
// Parameters:
//
// "old_trace" is a numeric value corresponding to the older sum of the squared diagonal elements.
//
// "new_trace" is a numeric value corresponding to the current sum of the squared diagonal elements.
//
// Return Value:
//
// returns a numeric value corresponding to the decimal difference between "old_trace" and "new_trace".

double PercentDifference(double old_trace, double new_trace)
{

    double percent_difference = (new_trace - old_trace) / old_trace;
    
    if (percent_difference < 0)
    {
    
        percent_difference = percent_difference * -1;
        
    }
    
    return percent_difference;
    
}

// "Electric_Field_x" calculates the x component of the electrostatic field using the
// finite difference formula.
//
// Parameters:
//
// "Nbox" is an integer value corresponding to half the total length of one side of the square
// encapsulating the parallel plates.
//
// "Nplate" is an integer value corresponding to the unit step size used to make the length
// of "Nbox" and "Nlplate".
//
// "Potential" is a pointer to an array holding the potential values from the final while loop iteration.
//
// Return Value:
//
// none

void Electric_Field_x(int Nbox, int Nplate, double* Potential)
{

    ofstream dataOut2("Electric_Field_x.txt");
    
    for (int i = -Nbox; i < Nbox; i++)
    {
    
        for (int j = -Nbox; j <= Nbox; j++)
        {
        
                dataOut2 << i << " " << j << " " << -(V(i + 1, j) - V(i, j)) << endl;
            
        }
        
    }
    
    dataOut2.close();
    
}

// "Electric_Field_y" calculates the y component of the electrostatic field using the
// finite difference formula.
//
// Parameters:
//
// "Nbox" is an integer value corresponding to half the total length of one side of the square
// encapsulating the parallel plates.
//
// "Nplate" is an integer value corresponding to the unit step size used to make the length
// of "Nbox" and "Nlplate".
//
// "Potential" is a pointer to an array holding the potential values from the final while loop iteration.
//
// Return Value:
//
// none

void Electric_Field_y(int Nbox, int Nplate, double* Potential)
{

    ofstream dataOut3("Electric_Field_y.txt");
    
    for (int i = -Nbox; i < Nbox; i++)
    {
    
        for (int j = -Nbox; j <= Nbox; j++)
        {
        
                dataOut3 << i << " " << j << " " << -(V(i, j + 1) - V(i, j)) << endl;
            
        }
        
    }
    
    dataOut3.close();
    
}

    
    