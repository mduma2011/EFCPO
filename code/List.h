#ifndef _LIST_H
#define _LIST_H

#include <assert.h>



template<class type>
class ArrayList
{	
	int size;
	type *vector;
	int last_index;

	void expand()
	{
		size = 2 * size;
		type *v = new type[size];
		for(int i = 0; i < last_index; i++) 
			v[i] = this->vector[i];
		delete vector;
		vector = v;
	}
	int BinarySearch(type data, int f, int l)
	{
		if(f == l)
		{
			if(data == vector[f]) return f;
			else if (data > vector[f]) return ~(f + 1);
			else return ~f;
		}
		if(f > l) 
		{
			if(last_index == 0) return ~f;
			if (data > vector[f]) 
				return ~(f + 1);
			return ~f;
		}
		int mid = (f + l) / 2;
		if(data == vector[mid]) return mid;
		else if(data < vector[mid]) return this->BinarySearch(data, f, mid - 1);
		else return this->BinarySearch(data, mid + 1, l);
	}
public:
	ArrayList()
	{
		last_index = 0;
		this->size = 16;
		vector = new type[size];
	}

	ArrayList(int initSize)
	{
		last_index = 0;
		this->size = initSize;
		vector = new type[size];
	}

	~ArrayList()
	{
		delete vector;
	}

	void Add(type input)
	{
		if(last_index >= size)
		{
			this->expand();
		}
		this->vector[last_index] = input;
		last_index++;
	}
	
	int BinarySearch(type data)
	{
		return BinarySearch(data, 0, last_index - 1);
	}

	int Count()
	{
		return last_index;
	}

	type Get(int index)
	{
		return vector[index];
	}

	void Set(int index, type value)
	{	
		vector[index] = value;
	}

	void Insert(int index, type input)
	{		
		if(last_index >= size)
		{
			this->expand();
		}
		for(int i = last_index; i > index; i--) vector[i] = vector[i - 1];
		vector[index] = input;
		last_index++;
	}

	void Delete(int index)
	{		
		for(int i = index; i < size - 1 ; i++) vector[i] = vector[i + 1];
		last_index--;
	}

	void DoClear()
	{		
		last_index = 0;
	}

	void DoClear(int count)
	{		
		last_index = count;
	}
};

#endif
