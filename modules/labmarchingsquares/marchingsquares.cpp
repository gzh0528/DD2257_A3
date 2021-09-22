/*********************************************************************
 *  Author  : Himangshu Saikia, Wiebke Koepp, Anke Friederici
 *  Init    : Monday, September 11, 2017 - 12:58:42
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labmarchingsquares/marchingsquares.h>
#include <inviwo/core/util/utilities.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo MarchingSquares::processorInfo_{
    "org.inviwo.MarchingSquares",  // Class identifier
    "Marching Squares",            // Display name
    "KTH Lab",                     // Category
    CodeState::Experimental,       // Code state
    Tags::None,                    // Tags
};
bool next_x(vec2 a,vec2 b)
{
    return a[0]<b[0];
}

const ProcessorInfo MarchingSquares::getProcessorInfo() const { return processorInfo_; }

MarchingSquares::MarchingSquares()
    : Processor()
    , inData("volumeIn")
    , meshIsoOut("meshIsoOut")
    , meshGridOut("meshGridOut")
    , propShowGrid("showGrid", "Show Grid")
    , propGridColor("gridColor", "Grid Lines Color", vec4(0.0f, 0.0f, 0.0f, 1.0f), vec4(0.0f),
                    vec4(1.0f), vec4(0.1f), InvalidationLevel::InvalidOutput,
                    PropertySemantics::Color)
    , propDeciderType("deciderType", "Decider Type")
    , propRandomSeed("seed", "Random Seed", 0, 0, std::mt19937::max())
    , propMultiple("multiple", "Iso Levels")
    , propIsoValue("isovalue", "Iso Value")
    , propIsoColor("isoColor", "Color", vec4(0.0f, 0.0f, 1.0f, 1.0f), vec4(0.0f), vec4(1.0f),
                   vec4(0.1f), InvalidationLevel::InvalidOutput, PropertySemantics::Color)
    , propNumContours("numContours", "Number of Contours", 1, 1, 50, 1)
    , propIsoTransferFunc("isoTransferFunc", "Colors", &inData) {
    // Register ports
    addPort(inData);
    addPort(meshIsoOut);
    addPort(meshGridOut);

    // Register properties
    addProperty(propShowGrid);
    addProperty(propGridColor);

    addProperty(propDeciderType);
    propDeciderType.addOption("asymptotic", "Asymptotic", 0);
    propDeciderType.addOption("random", "Random", 1);

    addProperty(propRandomSeed);
    propRandomSeed.setSemantics(PropertySemantics::Text);

    addProperty(propMultiple);

    propMultiple.addOption("single", "Single", 0);
    addProperty(propIsoValue);
    addProperty(propIsoColor);

    propMultiple.addOption("multiple", "Multiple", 1);
    addProperty(propNumContours);
    addProperty(propIsoTransferFunc);

    // The default transfer function has just two blue points
    propIsoTransferFunc.get().clear();
    propIsoTransferFunc.get().add(0.0f, vec4(0.0f, 0.0f, 1.0f, 1.0f));
    propIsoTransferFunc.get().add(1.0f, vec4(0.0f, 0.0f, 1.0f, 1.0f));
    propIsoTransferFunc.setCurrentStateAsDefault();

    util::hide(propGridColor, propRandomSeed, propNumContours, propIsoTransferFunc);

    propDeciderType.onChange([this]() {
        if (propDeciderType.get() == 1) {
            util::show(propRandomSeed);
        } else {
            util::hide(propRandomSeed);
        }
    });

    // Show the grid color property only if grid is actually displayed
    propShowGrid.onChange([this]() {
        if (propShowGrid.get()) {
            util::show(propGridColor);
        } else {
            util::hide(propGridColor);
        }
    });

    // Show options based on display of one or multiple iso contours
    propMultiple.onChange([this]() {
        if (propMultiple.get() == 0) {
            util::show(propIsoValue, propIsoColor);
            util::hide(propNumContours, propIsoTransferFunc);
        } else {
            util::hide(propIsoValue);
            util::show(propIsoColor, propNumContours);

            // TODO (Bonus): Comment out above if you are using the transfer function
            // and comment in below instead
            // util::hide(propIsoValue, propIsoColor);
            // util::show(propNumContours, propIsoTransferFunc);
        }
    });
}

void MarchingSquares::process() {
    if (!inData.hasData()) {
        return;
    }

    // Create a structured grid from the input volume
    auto vol = inData.getData();
    auto grid = ScalarField2::createFieldFromVolume(vol);

    // Extract the minimum and maximum value from the input data
    const double minValue = grid.getMinValue();
    const double maxValue = grid.getMaxValue();

    // Set the range for the isovalue to that minimum and maximum
    propIsoValue.setMinValue(minValue);
    propIsoValue.setMaxValue(maxValue);

    // You can print to the Inviwo console with Log-commands:
    LogProcessorInfo("This scalar field contains values between " << minValue << " and " << maxValue
                                                                  << ".");
    // You can also inform about errors and warnings:
    // LogProcessorWarn("I am warning about something"); // Will print warning message in yellow
    // LogProcessorError("I am letting you know about an error"); // Will print error message in red
    // (There is also LogNetwork...() and just Log...(), these display a different source,
    // LogProcessor...() for example displays the name of the processor in the workspace while
    // Log...() displays the identifier of the processor (thus with multiple processors of the
    // same kind you would not know which one the information is coming from

    // Get the definition of our structured grid with
    // - number of vertices in each dimension {nx, ny}
    const ivec2 nVertPerDim = grid.getNumVerticesPerDim();
    // - bounding box {xmin, ymin} - {xmax, ymax}
    const dvec2 bBoxMin = grid.getBBoxMin();
    const dvec2 bBoxMax = grid.getBBoxMax();
    // - cell size {dx, dy}
    const dvec2 cellSize = grid.getCellSize();

    // Values at the vertex positions can be accessed by the indices of the vertex
    // with index i ranging between [0, nx-1] and j in [0, ny-1]
    ivec2 ij = {0, 0};
    double valueAt00 = grid.getValueAtVertex(ij);
    LogProcessorInfo("The value at (0,0) is: " << valueAt00 << ".");

    // Initialize the output: mesh and vertices for the grid and bounding box
    auto gridmesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> gridvertices;

    auto indexBufferBBox = gridmesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    // bottomLeft to topLeft
    drawLineSegment(bBoxMin, vec2(bBoxMin[0], bBoxMax[1]), propGridColor.get(),
                    indexBufferBBox.get(), gridvertices);
    // topLeft to topRight
    drawLineSegment(vec2(bBoxMin[0], bBoxMax[1]), bBoxMax, propGridColor.get(),
                    indexBufferBBox.get(), gridvertices);
    // topRight to bottomRight
    drawLineSegment(bBoxMax, vec2(bBoxMax[0], bBoxMin[1]), propGridColor.get(),
                    indexBufferBBox.get(), gridvertices);
    // bottomRight to bottomLeft
    drawLineSegment(vec2(bBoxMax[0], bBoxMin[1]), bBoxMin, propGridColor.get(),
                    indexBufferBBox.get(), gridvertices);

    // Set the random seed to the one selected in the interface
    randGenerator.seed(static_cast<std::mt19937::result_type>(propRandomSeed.get()));
    // You can create a random sample between min and max with
    float minRand = 0.0;
    float maxRand = 1.0;
    float rand = randomValue(minRand, maxRand);
    LogProcessorInfo("The first random sample for seed " << propRandomSeed.get() << " between "
                                                         << minRand << " and " << maxRand << " is "
                                                         << rand << ".");

    // Properties are accessed with propertyName.get()
    if (propShowGrid.get()) {
        // TODO: Add grid lines of the given color
        //draw row lines
        auto indexBufferGridline = gridmesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
        for(int i=0;i<nVertPerDim[0]-1;i++)
        {
            for(int k=0;k<nVertPerDim[1];k++)
            {
            vec2 rv1=vec2(bBoxMin[0]+i*(cellSize[0]),bBoxMin[1]+k*(cellSize[1]));
            vec2 rv2=vec2(bBoxMin[0]+(i+1)*(cellSize[0]),bBoxMin[1]+(k)*(cellSize[1]));
            drawLineSegment(rv1, rv2, propGridColor.get(), indexBufferGridline.get(), gridvertices);
            }
        }
        for(int j=0;j<nVertPerDim[1]-1;j++)
        {
            for(int k=0;k<nVertPerDim[0];k++)
            {
                vec2 rv3=vec2(bBoxMin[0]+k*(cellSize[0]),bBoxMin[1]+j*(cellSize[1]));
                vec2 rv4=vec2(bBoxMin[0]+(k)*(cellSize[0]),bBoxMin[1]+(j+1)*(cellSize[1]));
                drawLineSegment(rv3, rv4, propGridColor.get(), indexBufferGridline.get(), gridvertices);
            }
            
        }
        // The function drawLineSegments creates two vertices at the specified positions,
        // that are placed into the Vertex vector defining our mesh.
        // An index buffer specifies which of those vertices should be grouped into to make up
        // lines/trianges/quads. Here two vertices make up a line segment.
        auto indexBufferGrid = gridmesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

        // Draw a line segment from v1 to v2 with a the given color for the grid
        vec2 v1 = vec2(0.5, 0.5);
        vec2 v2 = vec2(0.6, 0.7);
        drawLineSegment(v1, v2, propGridColor.get(), indexBufferGrid.get(), gridvertices);
    }

    // Set the created grid mesh as output
    gridmesh->addVertices(gridvertices);
    meshGridOut.setData(gridmesh);

    // TODO (Bonus) Gaussian filter
    // Our input is const (i.e. cannot be altered), but you need to compute smoothed data and write
    // it somewhere
    // Create an editable structured grid with ScalarField2 smoothedField =
    // ScalarField2(nVertPerDim, bBoxMin, bBoxMax - bBoxMin); Values can be set with
    // smoothedField.setValueAtVertex({0, 0}, 4.2);
    // and read again in the same way as before
    // smoothedField.getValueAtVertex(ij);
    // Initialize the output: mesh and vertices
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;
    std::cout<<propMultiple.get();
    if (propMultiple.get() == 0) {
        // TODO: Draw a single isoline at the specified isovalue (propIsoValue)
        auto indexBufferIsoline = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
        // and color it with the specified color (propIsoColor)
        for(int i=0;i<nVertPerDim[0]-1;i++)
        {
            for(int j=0;j<nVertPerDim[1]-1;j++)
            {
                vec2 v1=vec2(bBoxMin[0]+i*(cellSize[0]),bBoxMin[1]+j*(cellSize[1]));
                vec2 v3=vec2(bBoxMin[0]+i*(cellSize[0]),bBoxMin[1]+(j+1)*(cellSize[1]));
                vec2 v2=vec2(bBoxMin[0]+(i+1)*(cellSize[0]),bBoxMin[1]+j*(cellSize[1]));
                vec2 v4=vec2(bBoxMin[0]+(i+1)*(cellSize[0]),bBoxMin[1]+(j+1)*(cellSize[1]));
                
                vec2 a,b,c,d;
                //down
                std::vector<vec2> tp;
                std::cout<<v1[0]<<v1[1]<<std::endl;
                std::cout<<v2[0]<<v2[1]<<std::endl;
                float isoval=propIsoValue.get();
                if((grid.getValueAtVertex({i,j})-isoval)*((grid.getValueAtVertex({i+1,j})-isoval))<0)
                {
                    float y=v1[0]+(isoval-grid.getValueAtVertex({i,j}))*(v2[0]-v1[0])/(grid.getValueAtVertex({i+1,j})-grid.getValueAtVertex({i,j}));
                    a=vec2(y,v1[1]);
                    tp.push_back(a);
                }
                //left
                if((grid.getValueAtVertex({i,j})-isoval)*((grid.getValueAtVertex({i,j+1})-isoval))<0)
                {
                    float x=v1[1]+(isoval-grid.getValueAtVertex({i,j}))*(v3[1]-v1[1])/(grid.getValueAtVertex({i,j+1})-grid.getValueAtVertex({i,j}));
                    b=vec2(v1[0],x);
                    tp.push_back(b);
                }
                //right
                if((grid.getValueAtVertex({i+1,j})-isoval)*((grid.getValueAtVertex({i+1,j+1})-isoval))<0)
                {
                    float xx=v2[1]+(isoval-grid.getValueAtVertex({i+1,j}))*(v4[1]-v2[1])/(grid.getValueAtVertex({i+1,j+1})-grid.getValueAtVertex({i+1,j}));
                    c=vec2(v2[0],xx);
                    tp.push_back(c);
                }
                if((grid.getValueAtVertex({i,j+1})-isoval)*((grid.getValueAtVertex({i+1,j+1})-isoval))<0)
                {
                    float yy=v3[0]+(isoval-grid.getValueAtVertex({i,j+1}))*(v4[0]-v3[0])/(grid.getValueAtVertex({i+1,j+1})-grid.getValueAtVertex({i,j+1}));
                    d=vec2(yy,v3[1]);
                    tp.push_back(d);
                }
                if(tp.size()==2)
                {
                    std::cout<<"hahahha";
                    std::cout<<tp[0][0]<<tp[0][1];

                    drawLineSegment(tp[0], tp[1], propIsoColor.get(), indexBufferIsoline.get(), vertices);
                }
                else if (tp.size()==4)
                {
                    //“Asymptotic Decider”
                    std::sort(tp.begin(),tp.end(),next_x);
                    drawLineSegment(tp[0], tp[1], propIsoColor.get(), indexBufferIsoline.get(), vertices);
                    drawLineSegment(tp[2], tp[3], propIsoColor.get(), indexBufferIsoline.get(), vertices);
                    /*float fij=grid.getValueAtVertex({i,j});
                    float fij1=grid.getValueAtVertex({i,j+1});
                    float fi1j=grid.getValueAtVertex({i+1,j});
                    float fi1j1=grid.getValueAtVertex({i+1,j+1});
                    float fab=(fij*fi1j1-fi1j*fij1)/(fi1j1+fij-fij1-fi1j);
                    std::cout<<fab<<"heyI";
                    std::cout<<propIsoValue.get();
                    if(fab>=isoval)
                    {
                        drawLineSegment(tp[0], tp[1], propIsoColor.get(), indexBufferIsoline.get(), vertices);
                        drawLineSegment(tp[2], tp[3], propIsoColor.get(), indexBufferIsoline.get(), vertices);
                    }
                    else{
                        drawLineSegment(tp[0], tp[2], propIsoColor.get(), indexBufferIsoline.get(), vertices);
                        drawLineSegment(tp[1], tp[3], propIsoColor.get(), indexBufferIsoline.get(), vertices);
                    }*/
                }
                
            }
        }
    }

    else {
        // TODO: Draw the given number (propNumContours) of isolines between
        // the minimum and maximum value
        //First interpolate between minimum and maximum value
        //Bonus transfer function
        vec4 mincolor = propIsoTransferFunc.get().sample(0.0f);
        vec4 maxcolor = propIsoTransferFunc.get().sample(1.0f);
        std::vector<float> Nisovals;
        std::vector<vec4> NColors;
        std::cout<<mincolor[0]<<mincolor[1]<<mincolor[2]<<mincolor[3]<<std::endl;
        std::cout<<maxcolor[0]<<maxcolor[1]<<maxcolor[2]<<maxcolor[3]<<std::endl;
        for(int i=1;i<=propNumContours.get();i++)
        {
            Nisovals.push_back(minValue+(maxValue-minValue)*i/(propNumContours.get()+1));
            NColors.push_back(propIsoTransferFunc.get().sample(((maxValue-minValue)*i/(propNumContours.get()+1))/(maxValue-minValue)));
            std::cout<<((maxValue-minValue)*i/(propNumContours.get()+1))/(maxValue-minValue)<<std::endl;
        }
        
        for(int isn=0;isn<Nisovals.size();isn++)
        {
            std::cout<<NColors[isn][0]<<NColors[isn][1]<<NColors[isn][2]<<NColors[isn][3]<<std::endl;
            auto indexBufferIsoline = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
            // and color it with the specified color (propIsoColor)
            for(int i=0;i<nVertPerDim[0]-1;i++)
            {
                for(int j=0;j<nVertPerDim[1]-1;j++)
                {
                    vec2 v1=vec2(bBoxMin[0]+i*(cellSize[0]),bBoxMin[1]+j*(cellSize[1]));
                    vec2 v3=vec2(bBoxMin[0]+i*(cellSize[0]),bBoxMin[1]+(j+1)*(cellSize[1]));
                    vec2 v2=vec2(bBoxMin[0]+(i+1)*(cellSize[0]),bBoxMin[1]+j*(cellSize[1]));
                    vec2 v4=vec2(bBoxMin[0]+(i+1)*(cellSize[0]),bBoxMin[1]+(j+1)*(cellSize[1]));
                    
                    vec2 a,b,c,d;
                    //down
                    std::vector<vec2> tp;
                    float isoval=Nisovals[isn];
                    if((grid.getValueAtVertex({i,j})-isoval)*((grid.getValueAtVertex({i+1,j})-isoval))<0)
                    {
                        float y=v1[0]+(isoval-grid.getValueAtVertex({i,j}))*(v2[0]-v1[0])/(grid.getValueAtVertex({i+1,j})-grid.getValueAtVertex({i,j}));
                        a=vec2(y,v1[1]);
                        tp.push_back(a);
                    }
                    //left
                    if((grid.getValueAtVertex({i,j})-isoval)*((grid.getValueAtVertex({i,j+1})-isoval))<0)
                    {
                        float x=v1[1]+(isoval-grid.getValueAtVertex({i,j}))*(v3[1]-v1[1])/(grid.getValueAtVertex({i,j+1})-grid.getValueAtVertex({i,j}));
                        b=vec2(v1[0],x);
                        tp.push_back(b);
                    }
                    //right
                    if((grid.getValueAtVertex({i+1,j})-isoval)*((grid.getValueAtVertex({i+1,j+1})-isoval))<0)
                    {
                        float xx=v2[1]+(isoval-grid.getValueAtVertex({i+1,j}))*(v4[1]-v2[1])/(grid.getValueAtVertex({i+1,j+1})-grid.getValueAtVertex({i+1,j}));
                        c=vec2(v2[0],xx);
                        tp.push_back(c);
                    }
                    if((grid.getValueAtVertex({i,j+1})-isoval)*((grid.getValueAtVertex({i+1,j+1})-isoval))<0)
                    {
                        float yy=v3[0]+(isoval-grid.getValueAtVertex({i,j+1}))*(v4[0]-v3[0])/(grid.getValueAtVertex({i+1,j+1})-grid.getValueAtVertex({i,j+1}));
                        d=vec2(yy,v3[1]);
                        tp.push_back(d);
                    }
                    if(tp.size()==2)
                    {
                        //propIsoColor.get()
                        drawLineSegment(tp[0], tp[1], NColors[isn], indexBufferIsoline.get(), vertices);
                    }
                    else if (tp.size()==4)
                    {
                        //“Asymptotic Decider”
                        float fij=grid.getValueAtVertex({i,j});
                        float fij1=grid.getValueAtVertex({i,j+1});
                        float fi1j=grid.getValueAtVertex({i+1,j});
                        float fi1j1=grid.getValueAtVertex({i+1,j+1});
                        float fab=(fij*fi1j1-fi1j*fij1)/(fi1j1+fij-fij1-fi1j);
                        if(fab>=isoval)
                        {
                            drawLineSegment(tp[0], tp[1],  NColors[isn], indexBufferIsoline.get(), vertices);
                            drawLineSegment(tp[2], tp[3],  NColors[isn], indexBufferIsoline.get(), vertices);
                        }
                        else{
                            drawLineSegment(tp[0], tp[2],  NColors[isn], indexBufferIsoline.get(), vertices);
                            drawLineSegment(tp[1], tp[3],  NColors[isn], indexBufferIsoline.get(), vertices);
                        }
                    }
                    
                }
            }
        }
        // TODO (Bonus): Use the transfer function property to assign a color
        // The transfer function normalizes the input data and sampling colors
        // from the transfer function assumes normalized input, that means
        // vec4 color = propIsoTransferFunc.get().sample(0.0f);
        // is the color for the minimum value in the data
        // vec4 color = propIsoTransferFunc.get().sample(1.0f);
        // is the color for the maximum value in the data
    }

    // Note: It is possible to add multiple index buffers to the same mesh,
    // thus you could for example add one for the grid lines and one for
    // each isoline
    // Also, consider to write helper functions to avoid code duplication
    // e.g. for the computation of a single iso contour

    mesh->addVertices(vertices);
    meshIsoOut.setData(mesh);
}

float MarchingSquares::randomValue(const float min, const float max) const {
    return min + uniformReal(randGenerator) * (max - min);
}

void MarchingSquares::drawLineSegment(const vec2& v1, const vec2& v2, const vec4& color,
                                      IndexBufferRAM* indexBuffer,
                                      std::vector<BasicMesh::Vertex>& vertices) {
    // Add first vertex
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    // A vertex has a position, a normal, a texture coordinate and a color
    // we do not use normal or texture coordinate, but still have to specify them
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    // Add second vertex
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}

}  // namespace inviwo
