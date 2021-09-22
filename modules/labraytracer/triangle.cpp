/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Tuesday, October 17, 2017 - 10:24:30
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labraytracer/triangle.h>
#include <labraytracer/util.h>
#include <memory>

namespace inviwo {

Triangle::Triangle() {}

Triangle::Triangle(const vec3& v0, const vec3& v1, const vec3& v2, const vec3& uvw0,
                   const vec3& uvw1, const vec3& uvw2) {
    mVertices[0] = v0;
    mVertices[1] = v1;
    mVertices[2] = v2;
    mUVW[0] = uvw0;
    mUVW[1] = uvw1;
    mUVW[2] = uvw2;
}

bool point_in_triangle(const vec3& a, const vec3&  b, const vec3& c, const vec3& p)
{
    const vec3 pa = {a[0] - p[0], a[1] - p[1], a[2] - p[2]};
    const vec3 pb = { b[0] - p[0], b[1] - p[1], b[2] - p[2]};
    const vec3 pc = { c[0] - p[0], c[1] - p[1], c[2] - p[2] };

    const vec3 product_papb = cross(pa, pb);
    const vec3 product_pbpc = cross(pb, pc);
    const vec3 product_pcpa = cross(pc, pa);
    
    if (dot(product_papb, product_pbpc) > 0  &&  dot(product_papb, product_pcpa) > 0)
        return true;
    else
        return false;
}

bool Triangle::closestIntersection(const Ray& ray, double maxLambda,
                                   RayIntersection& intersection) const {
    // Programming TASK 1: Implement this method
    // Your code should compute the intersection between a ray and a triangle.
    //
    // If you detect an intersection, the return type should look similar to this:
    // if(rayIntersectsTriangle)
    //{
    //  intersection = RayIntersection(ray,shared_from_this(),lambda,ray.normal,uvw);
    //  return true;
    //}
    //
    // Hints :
    // Ray origin p_r : ray.getOrigin()
    // Ray direction t_r : ray.getDirection()
    // Compute the intersection point using ray.pointOnRay(lambda)
    const vec3 p_r=ray.getOrigin();
    const vec3 t_r =ray.getDirection();
    const vec3 t1=mVertices[1]-mVertices[0];
    const vec3 t2=mVertices[2]-mVertices[0];
    const vec3 n=cross(t1,t2);
    if (fabs(dot(t_r,n)) < Util::epsilon) {
        return false;
    }
    double lambda=dot(mVertices[0]-p_r,n)/(dot(t_r,n));
    if((lambda<0)||(lambda>maxLambda))
        return false;
    const vec3 p=ray.pointOnRay(lambda);
    // Inside triangle
    if(point_in_triangle(mVertices[0], mVertices[1], mVertices[2], p))
    {
        const vec3 uvw(0,0,0);
        Util::normalize(n);
        intersection = RayIntersection(ray,shared_from_this(),lambda,n,uvw);
        return  true;
    }
    return false;
}

bool Triangle::anyIntersection(const Ray& ray, double maxLambda) const {
    RayIntersection temp;
    return closestIntersection(ray, maxLambda, temp);
}

void Triangle::drawGeometry(std::shared_ptr<BasicMesh> mesh,
                            std::vector<BasicMesh::Vertex>& vertices) const {
    auto indexBuffer = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

    Util::drawLineSegment(mVertices[0], mVertices[1], vec4(0.2, 0.2, 0.2, 1), indexBuffer.get(),
                          vertices);
    Util::drawLineSegment(mVertices[1], mVertices[2], vec4(0.2, 0.2, 0.2, 1), indexBuffer.get(),
                          vertices);
    Util::drawLineSegment(mVertices[2], mVertices[0], vec4(0.2, 0.2, 0.2, 1), indexBuffer.get(),
                          vertices);
}

}  // namespace inviwo
