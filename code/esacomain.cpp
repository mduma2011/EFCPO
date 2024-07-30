/*
 * esacomain.cpp 
 * This class or file runs the ES ant colony optimisation algorithm
 */
#include <stdio.h>
#include <string.h>
#include "mex.h"
#include "Result.h"
#include "_OurACO.h"

//These methods are for analysis. 
void getArrayCharacteristics(const mxArray* array_instance);
void displaySubscript(const mxArray *array_ptr, mwSize index);
mxClassID analyseClass(const mxArray *array_ptr);

static void analyze_int32(const mxArray *array_ptr)
{
    mwSize total_num_of_elements, index;
    #if MX_HAS_INTERLEAVED_COMPLEX
        mxComplexInt32 *pc;
        mxInt32 *p;
        if(mxIsComplex(array_ptr)) 
        {
            pc = mxGetComplexInt32s(array_ptr);
            total_num_of_elements = mxGetNumberOfElements(array_ptr);

            for (index=0; index<total_num_of_elements; index++)  
            {
                mexPrintf("\t");
                displaySubscript(array_ptr, index);
                mexPrintf(" = %d + %di\n", (*pc).real,(*pc).imag);
                pc++;
            }
        }
        else 
        {
            p = mxGetInt32s(array_ptr);
            total_num_of_elements = mxGetNumberOfElements(array_ptr);

            for (index=0; index<total_num_of_elements; index++) 
            {
                mexPrintf("\t");
                displaySubscript(array_ptr, index);
                mexPrintf(" = %d\n", *p++);
            }
        }
    #else
        int  *pr, *pi;
        pr = (int *)mxGetData(array_ptr);
        pi = (int *)mxGetImagData(array_ptr);
        total_num_of_elements = mxGetNumberOfElements(array_ptr);

        for (index=0; index<total_num_of_elements; index++) 
        {
            mexPrintf("\t");
            displaySubscript(array_ptr, index);
            if (mxIsComplex(array_ptr)) 
            {
                mexPrintf(" = %d + %di\n", *pr++, *pi++);
            }
            else 
            {
                mexPrintf(" = %d\n", *pr++);
            }
        }
    #endif
}


static void analyze_int64(const mxArray *array_ptr)
{
    mwSize total_num_of_elements, index; 
    #if MX_HAS_INTERLEAVED_COMPLEX
        mxComplexInt64 *pc;
        mxInt64 *p;
        if(mxIsComplex(array_ptr)) 
        {
            pc = mxGetComplexInt64s(array_ptr);
            total_num_of_elements = mxGetNumberOfElements(array_ptr);

            for (index=0; index<total_num_of_elements; index++)  
            {
                mexPrintf("\t");
                displaySubscript(array_ptr, index);
                mexPrintf(" = %" FMT64 "d + %" FMT64 "di\n", (*pc).real,(*pc).imag);
                pc++;
            }
        }
        else 
        {
            p = mxGetInt64s(array_ptr);
            total_num_of_elements = mxGetNumberOfElements(array_ptr);

            for (index=0; index<total_num_of_elements; index++) 
            {
                mexPrintf("\t");
                displaySubscript(array_ptr, index);
                mexPrintf(" = %" FMT64 "d\n", *p++);
            }
        }
    #else
        int64_T  *pr, *pi;
        pr = (int64_T *)mxGetData(array_ptr);
        pi = (int64_T *)mxGetImagData(array_ptr);
        total_num_of_elements = mxGetNumberOfElements(array_ptr);

        for (index=0; index<total_num_of_elements; index++)  
        {
            mexPrintf("\t");
            displaySubscript(array_ptr, index);
            if (mxIsComplex(array_ptr)) 
            {
              mexPrintf(" = %" FMT64 "d + %" FMT64 "di\n", *pr++, *pi++);
            }
            else 
            {
              mexPrintf(" = %" FMT64 "d\n", *pr++);
            }
        }
    #endif
}

static void analyze_double(const mxArray *array_ptr)
{

    mwSize total_num_of_elements, index;
    #if MX_HAS_INTERLEAVED_COMPLEX
        mexPrintf("MX HAS INTERLEAVED COMPLEX\n");
        mxComplexDouble *pc;
        mxDouble *p;
        total_num_of_elements = mxGetNumberOfElements(array_ptr);
        if (mxIsComplex(array_ptr)) 
        {
            pc = mxGetComplexDoubles(array_ptr);
            for (index=0; index<total_num_of_elements; index++)  
            {
                mexPrintf("\t");
                displaySubscript(array_ptr, index);
                mexPrintf(" = %g + %gi\n",(*pc).real,(*pc).imag);
                pc++;
            }
        }
        else 
        {
            p = mxGetDoubles(array_ptr);
            for (index=0; index<total_num_of_elements; index++) 
            {
                mexPrintf("\t");
                displaySubscript(array_ptr, index);
                mexPrintf(" = %g\n", *p++);
            }
        }
    #else
        mexPrintf("NO INTERLEAVED COMPLEX\n");
        double *pr, *pi;
        total_num_of_elements = mxGetNumberOfElements(array_ptr);
        pr = mxGetPr(array_ptr);
        pi = mxGetPi(array_ptr);
        for (index = 0; index < total_num_of_elements; index++)  
        {
            mexPrintf("\t");
            displaySubscript(array_ptr, index);
            if (mxIsComplex(array_ptr)) 
            {
                mexPrintf(" = %g + %gi\n", *pr++, *pi++);
            }
            else 
            {
                mexPrintf(" = %g\n", *pr++);
            }
        }
    #endif
}

static void analyze_full(const mxArray *numeric_array_ptr)
{
    mxClassID  category;
    category = mxGetClassID(numeric_array_ptr);    
    switch (category)  
    {
        case mxINT32_CLASS:  
            analyze_int32(numeric_array_ptr);  
            break;      
        case mxINT64_CLASS:  
            analyze_int64(numeric_array_ptr);  
            break;      
        case mxDOUBLE_CLASS:
            analyze_double(numeric_array_ptr); 
            break;
        default: 
            break;
    }
}

/*
 * this method displays the row  and column. 
 */
void displaySubscript(const mxArray *array_ptr, mwSize index)
{
    mwSize     inner, subindex, total, d, q, number_of_dimensions;
    mwSize       *subscript;
    const mwSize *dims;

    number_of_dimensions = mxGetNumberOfDimensions(array_ptr);
    subscript = (mwSize *)mxCalloc(number_of_dimensions, sizeof(mwSize));
    dims = mxGetDimensions(array_ptr);

    mexPrintf("(");
    subindex = index;
    for (d = number_of_dimensions-1; ; d--) 
    { 
        for (total=1, inner=0; inner<d; inner++) 
        {
            total *= dims[inner];
        }
        subscript[d] = subindex / total;
        subindex = subindex % total;
        if (d == 0) {
            break;
        }
    }

    for (q = 0; q < number_of_dimensions-1; q++) 
    {
        mexPrintf("%d,", subscript[q] + 1);
    }
    mexPrintf("%d)", subscript[number_of_dimensions-1] + 1);

    mxFree(subscript);
}


/* Pass analyze_cell a pointer to a cell mxArray.  Each element
   in a cell mxArray is called a "cell"; each cell holds zero
   or one mxArray.  analyze_cell accesses each cell and displays
   information about it. */
static void analyze_cell(const mxArray *cell_array_ptr)
{
    mwSize total_num_of_cells;
    mwIndex index;
    const mxArray *cell_element_ptr;

    total_num_of_cells = mxGetNumberOfElements(cell_array_ptr);
    mexPrintf("total num of cells = %d\n", total_num_of_cells);
    mexPrintf("\n");

    /* Each cell mxArray contains m-by-n cells; Each of these cells
       is an mxArray. */
    for (index=0; index<total_num_of_cells; index++) 
    {
        mexPrintf("\n\n\t\tCell Element: ");
        displaySubscript(cell_array_ptr, index);
        mexPrintf("\n");
        cell_element_ptr = mxGetCell(cell_array_ptr, index);
        if (cell_element_ptr == NULL)
        {
            mexPrintf("\tEmpty Cell\n");
        }
        else 
        {
            /* Display a top banner. */
            mexPrintf("------------------------------------------------\n");
            getArrayCharacteristics(cell_element_ptr);
            analyseClass(cell_element_ptr);
            mexPrintf("\n");
        }
    }
	mexPrintf("\n");
}

/* Pass analyze_string a pointer to a char mxArray.  Each element
   in a char mxArray holds one 2-byte character (an mxChar); 
   analyze_string displays the contents of the input char mxArray
   one row at a time.  Since adjoining row elements are NOT stored in 
   successive indices, analyze_string has to do a bit of math to
   figure out where the next letter in a string is stored. */ 
static void analyze_string(const mxArray *string_array_ptr)
{
    char *buf;
    mwSize number_of_dimensions, buflen; 
    const mwSize *dims;
    mwSize d, page, total_number_of_pages, elements_per_page;

    /* Allocate enough memory to hold the converted string. */
    buflen = mxGetNumberOfElements(string_array_ptr) + 1;
    buf = (char*)mxCalloc(buflen, sizeof(char));

    /* Copy the string data from string_array_ptr and place it into buf. */
    if (mxGetString(string_array_ptr, buf, buflen) != 0) 
    {
        mexErrMsgIdAndTxt( "MATLAB:explore:invalidStringArray", "Could not convert string data.");
    }

    /* Get the shape of the input mxArray. */
    dims = mxGetDimensions(string_array_ptr);
    number_of_dimensions = mxGetNumberOfDimensions(string_array_ptr);

    elements_per_page = dims[0] * dims[1];
    /* total_number_of_pages = dims[2] x dims[3] x ... x dims[N-1] */
    total_number_of_pages = 1;
    for (d=2; d<number_of_dimensions; d++) 
    {
        total_number_of_pages *= dims[d];
    }

    for (page=0; page < total_number_of_pages; page++)
    {
        mwSize row;
        /* On each page, walk through each row. */
        for (row = 0; row < dims[0]; row++)  
        {
            mwSize column;
            mwSize index = (page * elements_per_page) + row;
            mexPrintf("\t");
            displaySubscript(string_array_ptr, index);
            mexPrintf(" ");

            /* Walk along each column in the current row. */
            for (column=0; column<dims[1]; column++) {
                mexPrintf("%c",buf[index]);
                index += dims[0];
            }
            mexPrintf("\n");
        }
    }
}

/* Determine the category (class) of the input array_ptr, and then
   branch to the appropriate analysis routine. */
mxClassID analyseClass(const mxArray *array_ptr)
{
    mxClassID  category;
    category = mxGetClassID(array_ptr);

    if (mxIsSparse(array_ptr)) {
       // analyze_sparse(array_ptr);
        mexPrintf("analyse sparse array\n");
    }
    else {
        switch (category) 
        {            
            case mxCHAR_CLASS:    
                 mexPrintf("char class\n");
                 analyze_string(array_ptr);   
                 break;
            //case mxLOGICAL_CLASS: analyze_logical(array_ptr);    break;
            //case mxSTRUCT_CLASS:  analyze_structure(array_ptr);  break;
            case mxCELL_CLASS: 
                 mexPrintf("cell class\n");
                 analyze_cell(array_ptr);       
                 break;
            case mxUNKNOWN_CLASS:
                mexWarnMsgIdAndTxt("MATLAB:explore:unknownClass", "Unknown class.");
                break;
            default:
                mexPrintf("full category class\n");
                analyze_full(array_ptr);
                break;
        }
    }
    return(category);
}

void getArrayCharacteristics(const mxArray* array_instance)
{
   
    //get the characteristics of array_instance
    const char* class_name;
    const mwSize* dims;
    char* shape_string;
    char* temp_string;
    mwSize c;
    mwSize number_of_dimensions;
    size_t length_of_shape;
      
    mexPrintf("------------------------------------------------------\n");
   
    //display the dimensions
    number_of_dimensions = mxGetNumberOfDimensions(array_instance);
    dims = mxGetDimensions(array_instance);
    
    shape_string = (char*)mxCalloc(3 * number_of_dimensions, sizeof(char));
    shape_string[0] = '\0';
    temp_string  = (char*)mxCalloc(64, sizeof(char));
       
    for(c = 0; c < number_of_dimensions; c++)
    {
      sprintf(temp_string, "%" FMT_SIZE_T "dx", dims[c]);
      strcat(shape_string, temp_string);
    }    
    mexPrintf("done all the copying to temp_string.\n");
    mexPrintf("let us see what it's copied: %s.\n", shape_string); //final outcome is something like this: 1x1x.
    
    length_of_shape = strlen(shape_string);
    
    //replace the last x with space
    shape_string[length_of_shape - 1] = '\0';    
    if(length_of_shape > 16)
    {
       sprintf(shape_string, "%" FMT_SIZE_T "u-D", number_of_dimensions); 
    }
    
    mexPrintf("Dimensions: %s\n", shape_string);
    
    
    //display the class category of mxArray
    class_name = mxGetClassName(array_instance);
                   
    mexPrintf("Class Name: %s%s\n", class_name, mxIsSparse(array_instance)? " (sparse)" : ""); 
    
    mexPrintf("------------------------------------------------------\n");
    mxFree(shape_string);
}

/**
 *  //display the characteristics of prhs[0],prhs[1],... or, prhs[n-1]
 *  getArrayCharacteristics(prhs[0]);    
 * // analyse the class
 * analyseClass(prhs[0]);
 **/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
    mexPrintf("nlhs: %d. nrhs: %d.\n", nlhs, nrhs);
    
    
    if (nrhs != 9) 
    {
        mexErrMsgIdAndTxt("MATLAB:esacomain:InvalidInput", "Invalid number of inputs to the esacomain MEX file. The number should be 9.");
    }
    
    if (nlhs != 4) 
    {
        mexErrMsgIdAndTxt("MATLAB:esacomain:InvalidOutput", "Invalid number of ouputs to the esacomain MEX file. The number should be 4.");
    }
    
      
    int number_of_ants = 0; 
    int number_of_iteration = 0; 
    double alpha = 0; 
    double rho = 0; 
    double beta = 0; 
    int max_Cl_count = 0; 
    int number_of_candidates_const = 0; 
    int number_of_candidates = 0;    
    char* fullpath = NULL;
    
    //read the input arrays
    double *pr;    
    for(int i = 0; i < nrhs; i++)
    {
        const mxArray* array_ptr = prhs[i];        
        if(mxIsDouble(array_ptr) || mxIsInt32(array_ptr) || mxIsInt64(array_ptr))
        {
            pr = mxGetPr(array_ptr);            
            switch(i)
            {
                case 0:                   
                    number_of_ants = *pr;                    
                    mexPrintf(" number_of_ants = %d : (*pr) : %g\n", number_of_ants);
                    break;
                case 1:                    
                    number_of_iteration = *pr;
                    mexPrintf(" number_of_iterations = %d\n", number_of_iteration);
                    break;   
                case 2:                    
                    alpha = (double)(*pr);
                    mexPrintf(" alpha = %0.3f\n", alpha);
                    break;   
                case 3:                    
                    rho = (double)(*pr);
                    mexPrintf(" rho = %0.3f\n", rho);
                    break;                    
                case 4:                    
                    beta = (double)(*pr);
                    mexPrintf(" beta = %0.3f\n", beta);
                    break;                    
                case 5:                    
                    max_Cl_count  = (int)(*pr);
                    mexPrintf(" max_Cl_count = %d\n", max_Cl_count);
                    break;                    
                case 6:                    
                    number_of_candidates_const = (int)(*pr);
                    mexPrintf(" number_of_candidates_const = %d\n", number_of_candidates_const);
                    break;
                case 7:                    
                    number_of_candidates = (int)(*pr);
                    mexPrintf(" number_of_candidates = %d\n", number_of_candidates);
                    break;
            }            
        }else if(i == nrhs-1)
        {
           //allocate memory to hold the converted string            
            mwSize buflength;
            buflength  =  mxGetNumberOfElements(array_ptr) + 1;
            fullpath = (char*)mxCalloc(buflength, sizeof(char));
            
            mxGetString(array_ptr, fullpath, buflength); 
            mexPrintf("Number of elems: %d .Full path : %s \n",buflength, fullpath);        
        }
        else
        {             
            mexErrMsgIdAndTxt("MATLAB:esacomain:InvalidInput", "Input values are not valid");
        }
    }
        
    Graph* graph = new Graph(fullpath);
    ESACO* esaco = new ESACO();
	
	Result res = esaco->MatlabRun(graph, number_of_ants, number_of_iteration, alpha, rho, beta, max_Cl_count, number_of_candidates_const, number_of_candidates);  
    int NumberOfNodes = graph->Dimension();
    mexPrintf("NumberOfNodes: %d.  Best Path Distance %d\n",NumberOfNodes, res.best_tour->Cost());
    
    double* antsMatrix;
    double* etaMatrix;
    double* besttour;
    double* besttourcost;
     
    for(int j = 0; j < nlhs; j++)
    {
       int index = 0;
       switch(j)
       {    
           case 0:               
               plhs[j] = mxCreateDoubleMatrix((mwSize)NumberOfNodes,(mwSize)number_of_ants, mxREAL);
               antsMatrix = mxGetPr(plhs[j]); 
               
               for(int r = 0; r < number_of_ants; r++)
               {                 
                  for(int c = 0; c < NumberOfNodes; c++)
                  {
                      antsMatrix[index++] = res.ants[r][c] + 1;                     
                  }                  
               }               
               break;
           case 1:               
               plhs[j] = mxCreateDoubleMatrix((mwSize)NumberOfNodes,(mwSize)NumberOfNodes, mxREAL);
               etaMatrix = mxGetPr(plhs[j]);
               
               for(int r = 0; r < NumberOfNodes; r++)
               {                 
                  for(int c = 0; c < NumberOfNodes; c++)
                  {                      
                      etaMatrix[index++] = (r == c)? 0 : res.ETA[r][c];                     
                  }                  
               }           
               break;                      
           case 2:               
                 plhs[j] = mxCreateDoubleMatrix((mwSize)NumberOfNodes, 1, mxREAL);  
                 besttour = mxGetPr(plhs[j]);
                 
                 for(int i =0; i < NumberOfNodes; i++)
                 {
                    besttour[i] = res.best_tour->Right(i) + 1;   
                 }                                       
              break;
            
           case 3:
                plhs[j] = mxCreateDoubleMatrix(1, 1, mxREAL);;
                besttourcost = mxGetPr(plhs[j]);
                besttourcost[0] = res.best_tour->Cost();
               break;
              
           default:
               break;
       }       
    }
    
    // Free some memory
    for(int i = 0; i < number_of_ants; i++)
    {
        delete res.ants[i];
    }  
    
    for(int i = 0; i < NumberOfNodes; i++)
    {
        delete res.ETA[i];
    }    
    delete res.ants;
    delete res.ETA;
    delete res.best_tour;
	delete graph;
    delete esaco; 
}
