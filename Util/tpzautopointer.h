//
// C++ Interface: tpzautopointer
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TPZAUTOPOINTER_H
#define TPZAUTOPOINTER_H

/**
This class implements a reference counter mechanism to administer a dynamically allocated object

@author Philippe R. B. Devloo
*/
template<class T>
class TPZAutoPointer{

template<class T2>
struct TPZReference
{
  T2 *fPointer;
  int fCounter;
  
  TPZReference()
  {
    fPointer = 0;
    fCounter = 1;
  }
  
  TPZReference(T2 *pointer)
  {
    fPointer = pointer;
    fCounter = 1;
  }
  
  ~TPZReference()
  {
    if(fPointer) delete fPointer;
  }
  
  
  /**
   * Increment the counter
   */
  void Increment()
  {
    fCounter++;
  }
  /**
   * Decrease the counter
   * If the counter is zero, delete myself
   */
  void Decrease()
  {
    fCounter--;
    if(fCounter <= 0) 
    {
      delete this;
    }
  }
  
};

/**
 * The object which contains the pointer and the reference count
 */
  TPZReference<T> *fRef;

public:
  /**
   * Creates an reference counted null pointer
   */
TPZAutoPointer()
{
  fRef = new TPZReference<T>();
}

/**
 * The destructor will delete the administered pointer if its reference count is zero
 */
~TPZAutoPointer()
{
  fRef->Decrease();
}
  
/**
 * This method will create an object which will administer the area pointed to by obj
 */
TPZAutoPointer(T *obj)
{
  fRef = new TPZReference<T>(obj);
}

/**
 * Share the pointer of the copy
 */
TPZAutoPointer(const TPZAutoPointer<T> &copy)
{
  fRef = copy.fRef;
  fRef->Increment();
}

TPZAutoPointer &operator=(const TPZAutoPointer<T> &copy)
{
  if(copy.fRef == fRef) return *this;
  copy.fRef->Increment();
  fRef->Decrease();
  fRef = copy.fRef;
  return *this;
}

T *operator->()
{
  return fRef->fPointer;
}

const T *operator->() const
{
  return fRef->fPointer;
}

};

#endif
