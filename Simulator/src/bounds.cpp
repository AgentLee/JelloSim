#include "bounds.h"

bool Bounds::Contains(Vector3f& pos)
{
    if( ( pos[0] > min[0] && pos[1] > min[1] && pos[2] > min[2] ) &&
        ( pos[0] < max[0] && pos[1] < max[1] && pos[2] < max[2] ) )
    {
        return true;
    }

    return false;
}

int Bounds::MaximumExtent() const
{
    Vector3f d = Diagonal();
    if (d[0] > d[1] && d[0] > d[2])
    {
        return 0;
    }
    else if (d[1] > d[2])
    {
        return 1;
    }
    else
    {
        return 2;
    }
}

Bounds Bounds::ApplyTransform(const Mat4f& T)
{
    //transform the min and max points of the bounding
    //box by the transformation tr
    Bounds b = Bounds();

    Point4f c1, c2, c3, c4, c5, c6, c7, c8;
    c1 << min[0], min[1], min[2], 1.0f;
    c2 << min[0], min[1], max[2], 1.0f;
    c3 << min[0], max[1], min[2], 1.0f;
    c4 << min[0], max[1], max[2], 1.0f;
    c5 << max[0], min[1], min[2], 1.0f;
    c6 << max[0], min[1], max[2], 1.0f;
    c7 << max[0], max[1], min[2], 1.0f;
    c8 << max[0], max[1], max[2], 1.0f;

    c1 = T * c1;
    c2 = T * c2;
    c3 = T * c3;
    c4 = T * c4;
    c5 = T * c5;
    c6 = T * c6;
    c7 = T * c7;
    c8 = T * c8;

    b = Union( b, c1 );
    b = Union( b, c2 );
    b = Union( b, c3 );
    b = Union( b, c4 );
    b = Union( b, c5 );
    b = Union( b, c6 );
    b = Union( b, c7 );
    b = Union( b, c8 );

    this->min = b.min;
    this->max = b.max;

    return b;
}

bool Bounds::Intersect(const Ray& r , float* t) const
{
    float t_n = -std::numeric_limits<T>::infinity();
    float t_f =  std::numeric_limits<T>::infinity();

    for(int i = 0; i < 3; i++)
    {
        //Ray parallel to slab check
        if(r.direction[i] == 0)
        {
            if(r.origin[i] < min[i] || r.origin[i] > max[i])
            {
                //The ray is parallel to the slab
                return false;
            }
        }

        //If not parallel, do slab intersect check
        float t0 = (min[i] - r.origin[i])/r.direction[i];
        float t1 = (max[i] - r.origin[i])/r.direction[i];

        if(t0 > t1)
        {
            float temp = t1;
            t1 = t0;
            t0 = temp;
        }

        if(t0 > t_n)
        {
            t_n = t0;
        }

        if(t1 < t_f)
        {
            t_f = t1;
        }
    }
    
    //making sure ray intersected the correct slab
    bool flag_in_box = false;

    if( max[0]-r.origin[0] > min[0] &&
        max[1]-r.origin[1] > min[1] &&
        max[2]-r.origin[2] > min[2] )
    {
        flag_in_box = true;
    }

    if(t_n < t_f)
    {
        float t = t_n > 0 ? t_n : t_f; //ensure ray moves in the positive direction
        if(t < 0 && !flag_in_box)
        {
            return false;
        }

        return true;
    }
    else //If t_near was greater than t_far, we did not hit the cube
    {
        return false;
    }

    return false;
}

Bounds Union(const Bounds& b1, const Bounds& b2)
{
    return Bounds( Point3f( std::min(b1.min[0], b2.min[0]),
                            std::min(b1.min[1], b2.min[1]),
                            std::min(b1.min[2], b2.min[2]) ),
                   Point3f( std::max(b1.max[0], b2.max[0]),
                            std::max(b1.max[1], b2.max[1]),
                            std::max(b1.max[2], b2.max[2]) ) );
}

Bounds Union(const Bounds& b1, const Point3f& p)
{
    return Bounds( Point3f( std::min(b1.min[0], p[0]),
                            std::min(b1.min[1], p[1]),
                            std::min(b1.min[2], p[2]) ),
                   Point3f( std::max(b1.max[0], p[0]),
                            std::max(b1.max[1], p[1]),
                            std::max(b1.max[2], p[2]) ) );
}

Bounds Union(const Bounds& b1, const Point4f& p)
{
    Point3f _p;
    _p << p[0], p[1], p[2];
    return Union( b1, _p );
}