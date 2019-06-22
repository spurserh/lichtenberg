
#ifndef COMMON_H
#define COMMON_H

#include <cmath>
#include <cassert>
#include <algorithm>

#if _WIN32
#pragma warning(disable:4250)	// Disable the "Inheritance via dominance" warning
#endif

#define PI						3.141592

typedef unsigned int	Size;
typedef unsigned Address;
	
template<typename Base>
struct Vec : public Base
{
	typedef typename Base::Type T;
	
	Vec()
	{
	}
	Vec(T x)
	{
		assert(this->mComponentCount==1);
		
		this->mComponents[0]=x;
	}
	Vec(T x,T y)
	{
		assert(this->mComponentCount==2);
		
		this->mComponents[0]=x;
		this->mComponents[1]=y;
	}
	Vec(T x,T y,T z)
	{
		assert(this->mComponentCount==3);
		
		this->mComponents[0]=x;
		this->mComponents[1]=y;
		this->mComponents[2]=z;
	}
	Vec(T x,T y,T z,T w)
	{
		assert(this->mComponentCount==4);
		
		this->mComponents[0]=x;
		this->mComponents[1]=y;
		this->mComponents[2]=z;
		this->mComponents[3]=w;
	}
	Vec(Vec const&o)
	{
		for(unsigned i=0;i<this->mComponentCount;++i)
			this->mComponents[i]=o.mComponents[i];
	}
#if 0
	Vec(VecBaseN<Base::Type,Base::mComponentCount> const&o)
	{
		for(unsigned i=0;i<this->mComponentCount;++i)
			this->mComponents[i]=o.mComponents[i];
	}
#endif
	
	Vec &operator=(Vec const&o)
	{
		for(unsigned i=0;i<this->mComponentCount;++i)
			this->mComponents[i]=o.mComponents[i];
		
		return *this;
	}
	
	bool operator==(Vec const&o)const
	{
		for(unsigned i=0;i<this->mComponentCount;++i)
		{
			if(this->mComponents[i]!=o.mComponents[i])
				return false;
		}
		return true;
	}
	
	bool operator!=(Vec const&o)const
	{
		return !((*this)==o);
	}
	
	bool operator<(Vec const&o)const
	{
		for(unsigned i=0;i<this->mComponentCount;++i)
		{
			if(this->mComponents[i]!=o.mComponents[i])
				return this->mComponents[i]<o.mComponents[i];
		}
		return false;
	}
	bool operator>(Vec const&o)const
	{
		for(unsigned i=0;i<this->mComponentCount;++i)
		{
			if(this->mComponents[i]!=o.mComponents[i])
				return this->mComponents[i]>o.mComponents[i];
		}
		return false;
	}
	bool operator<=(Vec const&o)const
	{
		if((*this)==o)
			return true;
		
		for(unsigned i=0;i<this->mComponentCount;++i)
		{
			if(this->mComponents[i]!=o.mComponents[i])
				return this->mComponents[i]<o.mComponents[i];
		}
		return false;
	}
	bool operator>=(Vec const&o)const
	{
		if((*this)==o)
			return true;
		
		for(unsigned i=0;i<this->mComponentCount;++i)
		{
			if(this->mComponents[i]!=o.mComponents[i])
				return this->mComponents[i]>o.mComponents[i];
		}
		return false;
	}
	
	Vec Cross(Vec const&o)const
	{
		assert((this->mComponentCount==3)&&(o.mComponentCount==3));
		return Vec(	this->y*o.z - this->z*o.y,
				   this->z*o.x - this->x*o.z,
				   this->x*o.y - this->y*o.x);
	}
	
	T Dot(Vec const&o)const
	{
		T accumulator=0;
		for(unsigned i=0;i<this->mComponentCount;++i)
			accumulator+=this->mComponents[i]*o.mComponents[i];
		return accumulator;
	}
	
	T &operator[](unsigned index)
	{
		assert(index<this->mComponentCount);
		return this->mComponents[index];
	}
	
	T const&operator[](unsigned index)const
	{
		assert(index<this->mComponentCount);
		return this->mComponents[index];
	}
	
	Vec operator-()const
	{
		Vec ret;
		for(unsigned i=0;i<this->mComponentCount;++i)
			ret[i]=-(*this)[i];
		return ret;
	}
	
	Vec &operator+=(Vec const&o)
	{
		for(unsigned i=0;i<this->mComponentCount;++i)
			this->mComponents[i]+=o.mComponents[i];
		
		return *this;
	}
	
	Vec &operator*=(float scalar)
	{
		for(unsigned i=0;i<this->mComponentCount;++i)
			this->mComponents[i]*=scalar;
		
		return *this;
	}
	
	Vec &operator-=(Vec const&o)
	{
		for(unsigned i=0;i<this->mComponentCount;++i)
			this->mComponents[i]-=o.mComponents[i];
		
		return *this;
	}
	
	Vec &operator/=(T scalar)
	{
		for(unsigned i=0;i<this->mComponentCount;++i)
			this->mComponents[i]/=scalar;
		
		return *this;
	}
	
	
	Vec operator+(Vec const&o)const
	{
		Vec ret;
		for(unsigned i=0;i<this->mComponentCount;++i)
			ret[i]=(*this)[i]+o[i];
		return ret;
	}
	
	Vec operator-(Vec const&o)const
	{
		Vec ret;
		for(unsigned i=0;i<this->mComponentCount;++i)
			ret[i]=(*this)[i]-o[i];
		return ret;
	}
	
	Vec operator*(typename Base::Type scalar)const
	{
		Vec ret;
		for(unsigned i=0;i<this->mComponentCount;++i)
			ret[i]=(*this)[i]*scalar;
		return ret;
	}
	
	Vec operator/(typename Base::Type scalar)const
	{
		Vec ret;
		for(unsigned i=0;i<this->mComponentCount;++i)
			ret[i]=(*this)[i]/scalar;
		return ret;
	}
	
	Vec operator*(Vec const&o)const
	{
		Vec ret;
		for(unsigned i=0;i<this->mComponentCount;++i)
			ret[i]=(*this)[i]*o[i];
		return ret;
	}
	
	Vec operator/(Vec const&o)const
	{
		Vec ret;
		for(unsigned i=0;i<this->mComponentCount;++i)
			ret[i]=(*this)[i]/o[i];
		return ret;
	}

	Vec &operator*=(Vec const&o)
	{
		for(unsigned i=0;i<this->mComponentCount;++i)
			(*this)[i]*=o[i];
		return (*this);
	}
	
	T Length()const
	{
		T accum=0;
		for(unsigned i=0;i<this->mComponentCount;++i)
			accum+=(*this)[i]*(*this)[i];
		return ::sqrt(accum);
	}
	
	Vec Normalized()const
	{
		return (*this)/Length();
	}
	
	T LargestComponent()const
	{
		assert(this->mComponentCount);
		T ret=this->mComponents[0];
		for(unsigned i=1;i<this->mComponentCount;++i)
		{
			if(this->mComponents[i]>ret)
				ret=this->mComponents[i];
		}
		return ret;
	}
    T InteriorAngleWith(Vec const&o)const
    {
        assert((this->mComponentCount==2)&&(o.mComponentCount==2));
        
        if((this->Length() == 0) || (o.Length() == 0))
            return 0.0f;
        
        Vec thisNormalized = this->Normalized();
        Vec oNormalized = o.Normalized();
        
        T thisAngle = ::atan2(oNormalized.y, oNormalized.x);
        T oAngle = ::atan2(thisNormalized.y, thisNormalized.x);
        
        float retAngle = ::fabs(thisAngle - oAngle);
		
        if(retAngle > PI)
            return (2*PI) - retAngle;
        
        return retAngle;
    }
};

template<typename T>
struct VecBase2
{
	typedef T Type;
	static const unsigned mComponentCount=2;
	
	VecBase2()
	{
		for(unsigned i=0;i<mComponentCount;++i)
			mComponents[i]=0;
	}
	
	T				&GetComponent(unsigned index)
	{
		return mComponents[index];
	}
	
	T		   const&GetComponentConst(unsigned index)const
	{
		return mComponents[index];
	}
	
	union
	{
		struct
		{
			T x,y;
		};
		struct
		{
			T width,height;
		};
		struct
		{
			T row,column;
		};
		T mComponents[2];
	};
};

template<typename T>
struct VecBase3
{
	typedef T Type;
	static const unsigned mComponentCount=3;
	
	VecBase3()
	{
		for(unsigned i=0;i<mComponentCount;++i)
			mComponents[i]=0;
	}
	
	T				&GetComponent(unsigned index)
	{
		return mComponents[index];
	}
	
	T		   const&GetComponentConst(unsigned index)const
	{
		return mComponents[index];
	}
	
	union
	{
		struct
		{
			T x,y,z;
		};
		struct
		{
			T width,height,depth;
		};
		struct
		{
			T row,column,plane;
		};
		struct
		{
			T red,green,blue;
		};
		T mComponents[3];
	};
};


template<typename T>
struct VecBase4
{
	typedef T Type;
	static const unsigned mComponentCount=4;
	
	VecBase4()
	{
		for(unsigned i=0;i<mComponentCount;++i)
			mComponents[i]=0;
	}
	
	T				&GetComponent(unsigned index)
	{
		return mComponents[index];
	}
	
	T		   const&GetComponentConst(unsigned index)const
	{
		return mComponents[index];
	}
	
	union
	{
		struct
		{
			T x,y,z,w;
		};
		struct
		{
			T red,green,blue,alpha;
		};
		T mComponents[4];
	};
};

template<typename T,unsigned N>
struct VecBaseN
{
	typedef T Type;
	static const unsigned mComponentCount=N;
	
	// Danger: Components are not set to 0
	// TODO: This may be fixable with template specialization (and without "memsetting")
	VecBaseN()
	{
	}
	
	T				&GetComponent(unsigned index)
	{
		return mComponents[index];
	}
	
	T		   const&GetComponentConst(unsigned index)const
	{
		return mComponents[index];
	}
	
	Vec<VecBase3<float> >	GetXYZ()
	{
		assert(mComponents>=3);
		return Vec<VecBase3<float> >(mComponents[0],mComponents[1],mComponents[2]);
	}
	
	T mComponents[N];
};


typedef Vec<VecBaseN<int,1> >		Vec1i;
typedef Vec<VecBase2<int> >		Vec2i;
typedef Vec<VecBase2<float> >		Vec2f;
typedef Vec<VecBase3<int> >		Vec3i;
typedef Vec<VecBase2<double> >	Vec2d;
typedef Vec<VecBase3<double> >	Vec3d;
typedef Vec<VecBase4<float> >	Vec4f;

template<typename T>
struct Quaternion : public Vec<VecBase4<T> >
{
	Quaternion(float x,float y,float z,float w)
		: Vec<VecBase4<T> >(x,y,z,w)
	{
	}

	Quaternion()
		: Vec<VecBase4<T> >(0.0f,0.0f,0.0f,1.0f)
	{
	}

	Quaternion(Quaternion const&o)
		: Vec<VecBase4<T> >(o)
	{
	}

	Vec<VecBase4<T> > &operator=(Vec<VecBase4<T> > const&o)
	{
		static_cast<Vec<VecBase4<T> > &>(*this) = static_cast<Vec<VecBase4<T> > const&>(o);
		
		return *this;
	}
	
	Quaternion operator *(Quaternion const&o)const
	{
		return Quaternion(
				((*this)[0] * o[3]) + ((*this)[3] * o[0]) + ((*this)[1] * o[2]) - ((*this)[2] * o[1]),
				((*this)[1] * o[3]) + ((*this)[3] * o[1]) + ((*this)[2] * o[0]) - ((*this)[0] * o[2]),
				((*this)[2] * o[3]) + ((*this)[3] * o[2]) + ((*this)[0] * o[1]) - ((*this)[1] * o[0]),
				((*this)[3] * o[3]) - ((*this)[0] * o[0]) - ((*this)[1] * o[1]) - ((*this)[2] * o[2])

			);
	}

	Quaternion operator *(Vec<VecBase3<T> > const&v)const
	{
		return Quaternion(
				  ((*this)[3] * v.x) + ((*this)[1] * v.z) - ((*this)[2] * v.y),
				  ((*this)[3] * v.y) + ((*this)[2] * v.x) - ((*this)[0] * v.z),
				  ((*this)[3] * v.z) + ((*this)[0] * v.y) - ((*this)[1] * v.x),
				- ((*this)[0] * v.x) - ((*this)[1] * v.y) - ((*this)[2] * v.z)
			);
	}

	Vec<VecBase3<T> > RotatePoint(Vec<VecBase3<T> > const&input)const
	{
	  Vec<VecBase3<T> > x;

	  Quaternion tmp, inv, f;

	  inv[0] = -(*this)[0]; 
	  inv[1] = -(*this)[1];
	  inv[2] = -(*this)[2]; 
	  inv[3] =  (*this)[3];

	  inv=inv.Normalized();

	  tmp = (*this) * input;
	  f = tmp * inv;

	  x.x = f[0];
	  x.y = f[1];
	  x.z = f[2];
	  
	  return x;
	}
};

typedef Quaternion<float> QuaternionF;
typedef Quaternion<float> QuaternionD;
	
struct Dimension
{
	inline Dimension():width(0),height(0){}
	inline Dimension(Size width,Size height):width(width),height(height){}
	inline Dimension(Dimension const&o):width(o.width),height(o.height){}
	inline Dimension(Vec2i const&o)
	: width(unsigned(o.width)),height(unsigned(o.height))
	{
		assert(o.width>=0);
		assert(o.height>=0);
	}
	
	inline Dimension&operator=(Dimension const&o)
	{
		width=o.width;
		height=o.height;
		return *this;
	}
	
	inline bool operator==(Dimension const&o)const
	{
		return (width==o.width)&&(height==o.height);
	}
	
	inline bool operator!=(Dimension const&o)const
	{
		return (width!=o.width)||(height!=o.height);
	}
	
	Size	width,height;
};
	
// This class should not be used within STL containers, ever!
template<typename T>
struct ArraySmartPtr{
	ArraySmartPtr():owner(false),ptr(0){}
	ArraySmartPtr(T *ptr)
	:	ptr(ptr),
	owner(ptr?true:false)
	{}	
	ArraySmartPtr(ArraySmartPtr const&o)
	:ptr(o.ptr),owner(o.owner)
	{
		o.owner=false;
	}
	~ArraySmartPtr()
	{
		if(owner)
		{
			delete []ptr;
			owner=false;
			ptr=0;
		}
	}
	ArraySmartPtr &operator=(T *ptr)
	{
		this->~ArraySmartPtr();
		this->ptr=ptr;
		owner=true;
		return *this;
	}
	ArraySmartPtr &operator=(ArraySmartPtr const&o)
	{
		this->~ArraySmartPtr();
		ptr=o.ptr;
		owner=o.owner;
		o.owner=false;
		return *this;
	}
	T *operator ->()const{return ptr;}
	operator T*()const{return ptr;}
	T *GetPtr()const{return ptr;}
	void DoNotDelete()const{owner=false;}
	
	T	*ptr;
private:
	mutable bool owner;
};

/**
 * SharedObjPtr can ONLY be used with heap-allocated objects deriving from SharedObject 
 *  (and thus allocated with its "new" function)
 *
 * PLEASE NOTE that there is currently no support for resolving cyclic references. 
 *
 * PLEASE NOTE that if SharedObjPtr is used to wrap the "this" pointer to an object during that object's
 *  constructor's execution, and that SharedObjPtr is then destroyed before the constructor exits, then 
 *  the object will be destroyed prematurely. 
 */
template<typename T>
struct SharedObjPtr
{
	SharedObjPtr()
	:	ptr(0)
	{
	}
	SharedObjPtr(T *ptr)
	:	ptr(ptr)
	{
		if(ptr)
			ptr->AddReference(true);
	}
	
	template<typename O>
	SharedObjPtr(SharedObjPtr<O> const&o)
	:	ptr(static_cast<T*>(o.GetPtr()))	// Try to convert the types (this enforces compatibility)
	{
		// We've just assigned our "ptr" using a dynamic cast, so if our ptr is 0 and theirs isn't, then
		//  that means that the types were incompatible.
		assert(!(o.GetPtr()&&(!ptr)));
		if(ptr)
			ptr->AddReference(false);
	}
	SharedObjPtr(SharedObjPtr const&o)
	:	ptr(o.ptr)
	{
		if(ptr)
			ptr->AddReference(false);
	}
	
	~SharedObjPtr()
	{
		if(ptr)
		{
			ptr->RemoveReference();
			ptr=0;
		}
	}
	
	template<typename O>
	SharedObjPtr &operator=(SharedObjPtr<O> const&o)
	{
		this->~SharedObjPtr();
		ptr=static_cast<T*>(o.GetPtr());	// Try to convert the types (this enforces compatibility)
		
		// We've just assigned our "ptr" using a dynamic cast, so if our ptr is 0 and theirs isn't, then
		//  that means that the types were incompatible.
		assert(!(o.GetPtr()&&(!ptr)));
		
		if(ptr)
			ptr->AddReference(false);
		return *this;
	}
	
	SharedObjPtr &operator=(SharedObjPtr const&o)
	{
		this->~SharedObjPtr();
		
		ptr=o.GetPtr();	// Try to convert the types (this enforces compatibility)
		
		// We've just assigned our "ptr" using a dynamic cast, so if our ptr is 0 and theirs isn't, then
		//  that means that the types were incompatible.
		assert(!(o.GetPtr()&&(!ptr)));
		
		if(ptr)
			ptr->AddReference(false);
		return *this;
	}
	
	bool operator<(SharedObjPtr const&o)const
	{
		Address thisA = reinterpret_cast<Address>(ptr);
		Address thisB = reinterpret_cast<Address>(o.ptr);

		return thisA<thisB;
	}
	
	T *operator ->()const
	{
		assert(ptr);
		return ptr;
	}
	operator T*()const
	{
		return ptr;
	}
	T *GetPtr()const
	{
		return ptr;
	}
	
private:
	T	*ptr;
};
	
class ISharedObject
{
	template<typename T>
	friend struct SharedObjPtr;
public:
	virtual
	~ISharedObject();

private:
	virtual void AddReference(bool fromRawPointer)const=0;
	virtual void RemoveReference()const=0;
};

struct SharedObject : public virtual ISharedObject
{
	template<typename T>
	friend struct SharedObjPtr;
	
	// Allocation of SharedObjects should be done with "new"
	SharedObject();
	
	/**
	 * No one except the garbage collector should ever delete a SharedObject. This is why the destructor is private. 
	 */
	virtual
	~SharedObject();
	
	
private:
		
	/**
	 * These are virtual because each module will end up with its own
	 *  object-tracker, and these need to refer to the correct one. 
	 *
	 * Also, each module may have its own allocator, and the same applies to that. 
	 */
	virtual void AddReference(bool fromRawPointer)const;
	virtual void RemoveReference()const;
	
	mutable unsigned													mReferenceCount;
};

template<typename Base>
struct Extrema
{
	Extrema()
	{
	}
	Extrema(Vec<Base> const&min,Vec<Base> const&max)
	: mMin(min), mMax(max)
	{
	}
	Extrema(Extrema const&o)
	: mMin(o.mMin), mMax(o.mMax)
	{
	}
	
	bool operator==(Extrema const&o)
	{
		for(unsigned i=0;i<Base::mComponentCount;++i)
		{
			if(mMin!=o.mMin)
				return false;
			if(mMax!=o.mMax)
				return false;
		}
		return true;
	}
	
	bool operator!=(Extrema const&o)
	{
		return !((*this)==o);
	}
	
	/**
	 * This function expands the Extrema as needed to enclose a point
	 */
	void DoEnclose(Vec<Base> const&pt)
	{
		for(unsigned i=0;i<Base::mComponentCount;++i)
		{
			if(pt.GetComponentConst(i)<mMin.GetComponentConst(i))
				mMin.GetComponent(i)=pt.GetComponentConst(i);
			
			if(pt.GetComponentConst(i)>mMax.GetComponentConst(i))
				mMax.GetComponent(i)=pt.GetComponentConst(i);
		}		
	}
	
	/**
	 * This function returns true if the given point is within the region,
	 *		false otherwise. 
	 */
	bool DoesEnclose(Vec<Base> const&pt)const
	{
		for(unsigned i=0;i<Base::mComponentCount;++i)
		{
			if(pt.GetComponentConst(i)<mMin.GetComponentConst(i))
				return false;
			
			if(pt.GetComponentConst(i)>mMax.GetComponentConst(i))
				return false;
		}
		return true;
	}
	
	/**
	 * This function finds the "union" of two sets of extrema, that is the Extrema that enclose both
	 */
	Extrema Union(Extrema const&o)const
	{
		Extrema ret;
		for(unsigned i=0;i<Base::mComponentCount;++i)
		{
			ret.mMin.GetComponent(i)=std::min(mMin.GetComponentConst(i), o.mMin.GetComponentConst(i));
			ret.mMax.GetComponent(i)=std::max(mMax.GetComponentConst(i), o.mMax.GetComponentConst(i));
		}
		return ret;
	}
	
	/**
	 * Gets the overlap (region covered by both extrema)
	 */
	Extrema GetOverlap(Extrema const&o)const
	{
		Extrema ret;
		for(unsigned i=0;i<Base::mComponentCount;++i)
		{
			ret.mMin.GetComponent(i)=(mMin.GetComponentConst(i)>o.mMin.GetComponentConst(i))?mMin.GetComponentConst(i):o.mMin.GetComponentConst(i);
			ret.mMax.GetComponent(i)=(mMax.GetComponentConst(i)<o.mMax.GetComponentConst(i))?mMax.GetComponentConst(i):o.mMax.GetComponentConst(i);
		}
		return ret;
	}

	bool DoesOverlap(Extrema const&o)const {
		const Extrema u = this->GetOverlap(o);
		for(unsigned i=0;i<Base::mComponentCount;++i)
		{
			if(u.mMin[i] >= u.mMax[i])
				return false;
		}
		return true;
	}
	
	/**
	 * Gets the size (width/height) of the region
	 */
	Vec<Base> GetSize()const
	{
		return mMax-mMin;
	}
	
	/**
	 * Returns true if the Extrema is valid (minimums <= maximums)
	 */
	bool IsValid()
	{
		return mMin <= mMax;
	}
	
	Vec<Base>	mMin,mMax;
};

typedef Extrema<VecBaseN<int, 1> >		Extrema1i;
typedef Extrema<VecBase3<int> >			Extrema3i;
typedef Extrema<VecBase2<int> >			Extrema2i;
typedef Extrema<VecBase2<float> >		Extrema2f;
typedef Extrema<VecBase2<double> >		Extrema2d;
typedef Extrema<VecBase3<double> >		Extrema3d;

class ConstantStringException : public std::exception
{
public:
	inline ConstantStringException(const char *what)
	: mWhat(what)
	{}
	
	inline const char *what()const throw()
	{
		return mWhat;
	}
	
private:
	
	const char *mWhat;
};


#endif//COMMON_H


