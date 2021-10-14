/*********************************************************************
 *  Author  : Anke Friederici
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 **********************************************************************/

#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/datastructures/volume/volumeram.h>
#include <labstreamlines/integrator.h>
#include <labutils/scalarvectorfield.h>
#include <labtopo/topology.h>
#include <labstreamlines/streamlineintegrator.h>
#include <labtopo/utils/gradients.h>

namespace inviwo {

const vec4 Topology::ColorsCP[6] = {
    vec4(1, 1, 0, 1),    // Saddle - Yellow
    vec4(1, 0, 0, 1),    // AttractingNode - Red
    vec4(0, 0, 1, 1),    // RepellingNode - Blue
    vec4(0.5, 0, 1, 1),  // AttractingFocus - Purple
    vec4(1, 0.5, 0, 1),  // RepellingFocus - Orange
    vec4(0, 1, 0, 1)     // Center - Green
};


// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo Topology::processorInfo_{
    "org.inviwo.Topology",    // Class identifier
    "Vector Field Topology",  // Display name
    "KTH Lab",                // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};

const ProcessorInfo Topology::getProcessorInfo() const { return processorInfo_; }

Topology::Topology()
    : Processor()
    , inData("inData")
    , outMesh("meshOut")
    , meshBBoxOut("meshBBoxOut")
    ,cellThreshold("Cell_Threshold", "ThresholdOfCellSize", {0.025,0.025}, {0.05,0.05})
    ,minThreshold("MinThreshold", "ThresholdofMin", 0.001,0.00001, 0.1,0.0001)
    , propStepSize("propStepSize", "Step size", 0.1, 0.001, 10)
    , propMaxSteps("propMaxSteps", "Max steps", 100, 0, 10000)
// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional), increment
// (optional)); propertyIdentifier cannot have spaces
{
    // Register Ports
    addPort(outMesh);
    addPort(inData);
    addPort(meshBBoxOut);

    // TODO: Register additional properties
    // addProperty(propertyName);
    addProperty(cellThreshold);
    addProperty(minThreshold);
    addProperty(propStepSize);
    addProperty(propMaxSteps);
}

void Topology::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);

    // Add a bounding box to the mesh
    const dvec2& BBoxMin = vectorField.getBBoxMin();
    const dvec2& BBoxMax = vectorField.getBBoxMax();
    auto bboxMesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> bboxVertices;
    auto indexBufferBBox = bboxMesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Bounding Box vertex 0
    vec4 black = vec4(0, 0, 0, 1);
    Integrator::drawNextPointInPolyline(BBoxMin, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMin[0], BBoxMax[1]), black, indexBufferBBox.get(),
                                        bboxVertices);
    Integrator::drawNextPointInPolyline(BBoxMax, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMax[0], BBoxMin[1]), black, indexBufferBBox.get(),
                                        bboxVertices);
    // Connect back to the first point, to make a full rectangle
    indexBufferBBox->add(static_cast<std::uint32_t>(0));
    bboxMesh->addVertices(bboxVertices);
    meshBBoxOut.setData(bboxMesh);

    // Initialize mesh, vertices and index buffers for seperatrices
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;
    // Either add all line segments to this index buffer (one large buffer, two consecutive points
    // make up one line), or use several index buffers with connectivity type strip.
    auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    // auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines,
    // ConnectivityType::Strip);

    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

    // TODO: Compute the topological skeleton of the input vector field.
    // Find the critical points and color them according to their type.
    // Integrate all separatrices.

    size2_t dims = vectorField.getNumVerticesPerDim();
    dvec2 cellsizeThres(0.005,0.005);
    std::cout<<vectorField.getBBoxMax().x<<"maxx"<<std::endl;
    std::cout<<vectorField.getBBoxMin().x<<"minx"<<std::endl;
    std::cout<<vectorField.getBBoxMax().y<<"maxy"<<std::endl;
    std::cout<<vectorField.getBBoxMin().y<<"miny"<<std::endl;
    // Looping through all values in the vector field.
    std::vector<dvec2> res;
    for (size_t j = 0; j < dims[1]; ++j) {
        for (size_t i = 0; i < dims[0]; ++i) {
            //change of sign test
            //vector<dvec2> res=extractCriticalPoints(vectorField, dims,dvec2(i,j));
            if((i+1<dims[0])&&(j+1<dims[1]))
            {
                std::vector<dvec2> gridPos;
                std::vector<dvec2> gridres;
                gridPos.push_back(vectorField.getPositionAtVertex(size2_t(i,j)));
                gridPos.push_back(vectorField.getPositionAtVertex(size2_t(i,j+1)));
                gridPos.push_back(vectorField.getPositionAtVertex(size2_t(i+1,j)));
                gridPos.push_back(vectorField.getPositionAtVertex(size2_t(i+1,j+1)));
                extractCriticalPoints(vectorField, gridPos,res,cellsizeThres);
            }
            
            
        }
    }
    for(auto cp:res)
    {
        int color_=criticalPointType(vectorField, res, cp);
        if(color_!=-1)
        {
            vec4 col_=ColorsCP[color_];
            Integrator::drawPoint(cp, col_, indexBufferPoints.get(), vertices);
        }
        else
        {
            Integrator::drawPoint(cp, black, indexBufferPoints.get(), vertices);
        }
        
    }
    for(auto cp:res)
    {
        //Topology::separatrice(vectorField, indexBufferPoints,mesh,vertices,cp);
    }
    
    // Other helpful functions
    // dvec2 pos = vectorField.getPositionAtVertex(size2_t(i, j));
    // Computing the jacobian at a position
    // dmat2 jacobian = vectorField.derive(pos);
    // Doing the eigen analysis
    // auto eigenResult = util::eigenAnalysis(jacobian);
    // The result of the eigen analysis has attributed eigenvaluesRe eigenvaluesIm and
    // eigenvectors
    
    // Accessing the colors
    vec4 colorCenter = ColorsCP[static_cast<int>(TypeCP::Center)];
    
    mesh->addVertices(vertices);
    outMesh.setData(mesh);
}

void Topology::drawLineSegment(const dvec2& v1, const dvec2& v2, const vec4& color,
                               IndexBufferRAM* indexBuffer,
                               std::vector<BasicMesh::Vertex>& vertices) {
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}

void Topology::extractCriticalPoints(const VectorField2& vectorField,std::vector<dvec2> gridPos, std::vector<dvec2>& res,dvec2 cellT)
{
    
    if(gridPos.size()!=4)
        return;
    /*if((!vectorField.isInside(gridPos[0]))||(!vectorField.isInside(gridPos[1]))||(!vectorField.isInside(gridPos[2]))||(!vectorField.isInside(gridPos[3])))
    {
        std::cout<<"not in side"<<gridPos[0]<<gridPos[1]<<gridPos[2]<<gridPos[3]<<std::endl;
        return;
    }*/
    
    /*
     If we recursively achieve a very small grid, then we collect results and
     exit the recursion*/
    if((gridPos[1].y-gridPos[0].y<=cellT.y)&&(gridPos[2].x-gridPos[0].x <= cellT.x))
    {
        //change of sign test begins
        dvec2 vij = vectorField.interpolate(gridPos[0]);
        dvec2 vij1 = vectorField.interpolate(gridPos[1]);
        dvec2 vi1j = vectorField.interpolate(gridPos[2]);
        dvec2 vi1j1 = vectorField.interpolate(gridPos[3]);
        int x1=vij.x>0?1:0;
        int x2=vij1.x>0?1:0;
        int x3=vi1j.x>0?1:0;
        int x4=vi1j1.x>0?1:0;
        int y1=vij.y>0?1:0;
        int y2=vij1.y>0?1:0;
        int y3=vi1j.y>0?1:0;
        int y4=vi1j1.y>0?1:0;
        if(((x1==x2)&&(x2==x3)&&(x3==x4))||((y1==y2)&&(y2==y3)&&(y3==y4)))
        {
            return;
        }
        /*'
         1. linear interpolation four edges
         2. v.x==0 then ->vx
         3. v.y==0 then ->vy
         4. test if the zero point is on the boundary
         5. get intersectionpoints of two lines and return
         */
        std::vector<dvec2> vx;
        std::vector<dvec2> vy;
        //bottom
        if(vij.x*vij1.x<0)
        {
            dvec2 px;
            px=dvec2(gridPos[0].x,gridPos[0].y+(0-vij.x)*(gridPos[1].y-gridPos[0].y)/(vij1.x-vij.x));
            if((vectorField.interpolate(px).x<minThreshold.get())&&(vectorField.interpolate(px).x>-minThreshold.get()))
            {
                LogProcessorInfo("px"<<vectorField.interpolate(px).x);
                vx.push_back(px);
            }
        }
        //up
        if(vi1j.x*vi1j1.x<0)
        {
            dvec2 px;
            px=dvec2(gridPos[2].x,gridPos[2].y+(0-vi1j.x)*(gridPos[3].y-gridPos[2].y)/(vi1j1.x-vi1j.x));
            if((vectorField.interpolate(px).x<minThreshold.get())&&(vectorField.interpolate(px).x>-minThreshold.get()))
            {
                LogProcessorInfo("px"<<vectorField.interpolate(px).x);
                vx.push_back(px);
            }
        }
        //left
        if(vij.x*vi1j.x<0)
        {
            dvec2 px;
            px=dvec2(gridPos[0].x+(0-vij.x)*(gridPos[2].x-gridPos[0].x)/(vi1j.x-vij.x),gridPos[0].y);
            if((vectorField.interpolate(px).x<minThreshold.get())&&(vectorField.interpolate(px).x>-minThreshold.get()))
            {
                LogProcessorInfo("px"<<vectorField.interpolate(px).x);
                vx.push_back(px);
            }
        }
        //right
        if(vij1.x*vi1j1.x<0)
        {
            dvec2 px;
            px=dvec2(gridPos[1].x+(0-vij1.x)*(gridPos[3].x-gridPos[1].x)/(vi1j1.x-vij1.x),gridPos[1].y);
            if((vectorField.interpolate(px).x<minThreshold.get())&&(vectorField.interpolate(px).x>-minThreshold.get()))
            {
                LogProcessorInfo("px"<<vectorField.interpolate(px).x);
                vx.push_back(px);
            }
        }
        
        //bottom
        if(vij.y*vij1.y<0)
        {
            dvec2 py;
            py=dvec2(gridPos[0].x,gridPos[0].y+(0-vij.y)*(gridPos[1].y-gridPos[0].y)/(vij1.y-vij.y));
            if((vectorField.interpolate(py).y<minThreshold.get())&&(vectorField.interpolate(py).y>-minThreshold.get()))
            {
                vy.push_back(py);
            }
        }
        //up
        if(vi1j.y*vi1j1.y<0)
        {
            dvec2 py;
            py=dvec2(gridPos[2].x,gridPos[2].y+(0-vi1j.y)*(gridPos[3].y-gridPos[2].y)/(vi1j1.y-vi1j.y));
            if((vectorField.interpolate(py).y<minThreshold.get())&&(vectorField.interpolate(py).y>-minThreshold.get()))
            {
                vy.push_back(py);
            }
        }
        //left
        if(vij.y*vi1j.y<0)
        {
            dvec2 py;
            py=dvec2(gridPos[0].x+(0-vij.y)*(gridPos[2].x-gridPos[0].x)/(vi1j.y-vij.y),gridPos[0].y);
            if((vectorField.interpolate(py).y<minThreshold.get())&&(vectorField.interpolate(py).y>-minThreshold.get()))
            {
                vy.push_back(py);
            }
        }
    
        //right
        if(vij1.y*vi1j1.y<0)
        {
            dvec2 py;
            py=dvec2(gridPos[1].x+(0-vij1.y)*(gridPos[3].x-gridPos[1].x)/(vi1j1.y-vij1.y),gridPos[1].y);
            if((vectorField.interpolate(py).y<minThreshold.get())&&(vectorField.interpolate(py).y>-minThreshold.get()))
            {
                vy.push_back(py);
            }
        }
        //hit the bound
        for(int ii=0;ii<vx.size();ii++)
        {
            dvec2 vii=vectorField.interpolate(vx[ii]);
            dmat2 jaco=vectorField.derive(vii);
            if((glm::determinant(jaco))&&(vii.x<minThreshold.get())&&(vii.x>-minThreshold.get())&&(vii.y<minThreshold.get())&&(vii.y>-minThreshold.get()))
            {
                res.push_back(vx[ii]);
                vx.erase(vx.begin()+ii);
            }
        }
        for(int jj=0;jj<vy.size();jj++)
        {
            dvec2 vjj=vectorField.interpolate(vy[jj]);
            dmat2 jaco=vectorField.derive(vjj);
            if((glm::determinant(jaco))&&(vjj.x<minThreshold.get())&&(vjj.x>-minThreshold.get())&&(vjj.y<minThreshold.get())&&(vjj.y>-minThreshold.get()))
            {
                res.push_back(vy[jj]);
                vy.erase(vy.begin()+jj);
            }
        }
        /*
         This is where we return the intersection point
         */
        if((vx.size()==2)&&(vy.size()==2))
        {
            dvec2 npos=Topology::lineLineIntersection(vx[0], vx[1], vy[0], vy[1]);
            if(npos.x!=FLT_MAX)
            {
                dmat2 jacobian = vectorField.derive(npos);
                //compute determinate
                float det=glm::determinant(jacobian);//jacobian[0][0]*jacobian[1][1]-jacobian[0][1]*jacobian[1][0];
                if(det!=0.0)
                    res.push_back(npos);
            }
        }
        return;
    }
    
    
    /*change of sign test*/
    dvec2 vij = vectorField.interpolate(gridPos[0]);
    dvec2 vij1 = vectorField.interpolate(gridPos[1]);
    dvec2 vi1j = vectorField.interpolate(gridPos[2]);
    dvec2 vi1j1 = vectorField.interpolate(gridPos[3]);
    int x1=vij.x>0?1:0;
    int x2=vij1.x>0?1:0;
    int x3=vi1j.x>0?1:0;
    int x4=vi1j1.x>0?1:0;
    int y1=vij.y>0?1:0;
    int y2=vij1.y>0?1:0;
    int y3=vi1j.y>0?1:0;
    int y4=vi1j1.y>0?1:0;
    if(((x1==x2)&&(x2==x3)&&(x3==x4))||((y1==y2)&&(y2==y3)&&(y3==y4)))
    {
        return;
    }
    
    /*
     Recursively decomposing the big grid
     */
    std::vector<dvec2> DecompPos(5,dvec2(0));
    DecompPos[0]=(gridPos[0]+gridPos[1])/2.0;//bottom
    DecompPos[1]=(gridPos[0]+gridPos[2])/2.0;//left
    DecompPos[2]=(gridPos[1]+gridPos[3])/2.0;//right
    DecompPos[3]=(gridPos[2]+gridPos[3])/2.0;//up
    DecompPos[4]=(gridPos[0]+gridPos[3])/2.0;//center
    //up left
    std::vector<dvec2>upleft={DecompPos[1],DecompPos[4],gridPos[2],DecompPos[3]};
    extractCriticalPoints(vectorField,upleft,res,cellT);
 
    //bottom left
    std::vector<dvec2>bottomleft={gridPos[0],DecompPos[0],DecompPos[1],DecompPos[4]};
    extractCriticalPoints(vectorField,bottomleft,res,cellT);
 
    //up right
    std::vector<dvec2>upright={DecompPos[4],DecompPos[2],DecompPos[3],gridPos[3]};
    extractCriticalPoints(vectorField,upright,res,cellT);

    //bottom right
    std::vector<dvec2>bottomright={DecompPos[0],gridPos[1],DecompPos[4],DecompPos[2]};
    extractCriticalPoints(vectorField,bottomright,res,cellT);


   
}

//not working
void Topology::separatrice(const VectorField2& vectorField,
                 std::shared_ptr<IndexBufferRAM> indexBufferPoints,
                 std::shared_ptr<BasicMesh> mesh,
                 std::vector<BasicMesh::Vertex>& vertices,dvec2 pos)
{
   
    dmat2 jacobian = vectorField.derive(pos);
    auto eigenResult = util::eigenAnalysis(jacobian);
    float MinT=minThreshold.get();
    dvec2 np1(pos.x+MinT*eigenResult.eigenvectors[0].x,pos.y+MinT*eigenResult.eigenvectors[0].y);
    dvec2 np2(pos.x-MinT*eigenResult.eigenvectors[0].x,pos.y-MinT*eigenResult.eigenvectors[0].y);
    dvec2 np3(pos.x+MinT*eigenResult.eigenvectors[1].x,pos.y+MinT*eigenResult.eigenvectors[1].y);
    dvec2 np4(pos.x-MinT*eigenResult.eigenvectors[1].x,pos.y-MinT*eigenResult.eigenvectors[1].y);
    Topology::drawStreamline(vectorField, np1, indexBufferPoints, mesh, vertices,propStepSize.get());
    Topology::drawStreamline(vectorField, np1, indexBufferPoints, mesh, vertices,-propStepSize.get());
    Topology::drawStreamline(vectorField, np1, indexBufferPoints, mesh, vertices,propStepSize.get());
    Topology::drawStreamline(vectorField, np1, indexBufferPoints, mesh, vertices,-propStepSize.get());
    
}
int Topology::drawStreamline(
        const VectorField2& vectorField,
        dvec2 startPoint,
        std::shared_ptr<IndexBufferRAM> indexBufferPoints,
        std::shared_ptr<BasicMesh> mesh,
                   std::vector<BasicMesh::Vertex>& vertices,double stepsize)
{
    // Draw start point
    // TODO: Create one stream line from the given start point
    auto indexBufferLine = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    vec4 color(0, 0, 1, 1);
    Integrator::drawNextPointInPolyline(startPoint, color, indexBufferLine.get(), vertices);
    dvec2 x0 = startPoint;
    int stepsTaken = 0;
    for (;;) {
        dvec2 x1 = Integrator::RK4_norm(vectorField, x0, stepsize);
        if (!vectorField.isInside(x1)) {
            LogProcessorInfo("Stop: domain");
            break;
        }

        if (vectorField.interpolate(x1) == dvec2(0)) {
            LogProcessorInfo("Stop: zero");
            break;
        }
        Integrator::drawNextPointInPolyline(x1, color, indexBufferLine.get(), vertices);
        x0 = x1;
        stepsTaken++;
    }

    return stepsTaken;
}

int Topology::criticalPointType(const VectorField2 &vectorField,std::vector<dvec2>& res,dvec2 pos)
{
    
            /*
             (-minThreshold,minThreshold) -> zero
             */
            dmat2 jacobian = vectorField.derive(pos);
            auto eigenResult = util::eigenAnalysis(jacobian);
        
            //compute determinate
            float det=glm::determinant(jacobian);
            if((eigenResult.eigenvaluesRe[0]*eigenResult.eigenvaluesRe[1]<0.0)&&((eigenResult.eigenvaluesIm[0]<minThreshold.get()&&(eigenResult.eigenvaluesIm[0]>-minThreshold.get()))&&((eigenResult.eigenvaluesIm[1]<minThreshold.get()&&(eigenResult.eigenvaluesIm[1]>-minThreshold.get())))))
            {
                
                LogProcessorInfo("Saddle Point")
                return 0;
            }
            if((eigenResult.eigenvaluesRe[0]>0)&& (eigenResult.eigenvaluesRe[1]>0)&&((eigenResult.eigenvaluesIm[0]<minThreshold.get()&&(eigenResult.eigenvaluesIm[0]>-minThreshold.get()))&&((eigenResult.eigenvaluesIm[1]<minThreshold.get()&&(eigenResult.eigenvaluesIm[1]>-minThreshold.get())))))
            {
                LogProcessorInfo("Repelling Node")
                return 2;
            }
    if((eigenResult.eigenvaluesRe[0]<0)&& (eigenResult.eigenvaluesRe[1]<0)&&((eigenResult.eigenvaluesIm[0]<minThreshold.get()&&(eigenResult.eigenvaluesIm[0]>-minThreshold.get()))&&((eigenResult.eigenvaluesIm[1]<minThreshold.get()&&(eigenResult.eigenvaluesIm[1]>-minThreshold.get())))))
            {
                LogProcessorInfo("Attracting Node")
                return 1;
            }
            if ((eigenResult.eigenvaluesRe[0]==eigenResult.eigenvaluesRe[1])&& ((eigenResult.eigenvaluesRe[1]<minThreshold.get())&&(eigenResult.eigenvaluesRe[1]>-minThreshold.get()))&&((eigenResult.eigenvaluesIm[0]==-eigenResult.eigenvaluesIm[1])&&(eigenResult.eigenvaluesIm[1]!=0.0))) {
                LogProcessorInfo("Center")
                return 5;
            }
            
            if ((eigenResult.eigenvaluesRe[0]==eigenResult.eigenvaluesRe[1])&& (eigenResult.eigenvaluesRe[1]<0.0)&&((eigenResult.eigenvaluesIm[0]==-eigenResult.eigenvaluesIm[1])&&(eigenResult.eigenvaluesIm[1]!=0.0))) {
                LogProcessorInfo("Attracting Focus")
                return 3;
            }
            if ((eigenResult.eigenvaluesRe[0]==eigenResult.eigenvaluesRe[1])&& (eigenResult.eigenvaluesRe[1]>0.0)&&((eigenResult.eigenvaluesIm[0]==-eigenResult.eigenvaluesIm[1])&&(eigenResult.eigenvaluesIm[1]!=0.0))) {
                LogProcessorInfo("Repelling Focus")
                return 4;
            }
    return -1;
}
dvec2 Topology::lineLineIntersection(dvec2 A, dvec2 B, dvec2 C, dvec2 D)
{
    // Line AB
    double a1 = B.y - A.y;
    double b1 = A.x - B.x;
    double c1 = a1*(A.x) + b1*(A.y);
    // Line CD
    double a2 = D.y - C.y;
    double b2 = C.x - D.x;
    double c2 = a2*(C.x)+ b2*(C.y);
  
    double determinant = a1*b2 - a2*b1;
  
    if (determinant == 0)
    {
        return dvec2(FLT_MAX, FLT_MAX);
    }
    else
    {
        double x = (b2*c1 - b1*c2)/determinant;
        double y = (a1*c2 - a2*c1)/determinant;
        return dvec2(x, y);
    }
}

}  // namespace inviwo
