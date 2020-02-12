#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color = new pixel[width*height];
    state.image_depth=new float[width*height];
    std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;

    // memset(state.image_color, -1, width*height*sizeof(pixel));
    for(int i = 0; i < height*width; i++)
    {
      state.image_color[i] = make_pixel(0,0,0);
      state.image_depth[i] = 1;
    }
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    // std::cout<<"TODO: implement rendering."<<std::endl;
    // int halfHeight = state.image_height/2;
    // for(int i = 0; i < state.image_width; i++)
    // {
    //   state.image_color[halfHeight*state.image_width+i] = make_pixel(0,255,255);
    // }
    // printf("num vertices %d\n", state.num_vertices);
    // printf("floats per vertex %d\n", state.floats_per_vertex);
    // printf("num triangles %d\n", state.num_triangles);
    struct data_vertex tempVertex;
    struct data_geometry tempGeometry;
    tempGeometry.data = new float[MAX_FLOATS_PER_VERTEX];

    struct data_fragment tempFragment;
    struct data_output tempOutput;
    tempFragment.data = new float[MAX_FLOATS_PER_VERTEX];



    int numTriangles = 0;
    if(type == render_type::triangle)
    {
      numTriangles = state.num_vertices/3;
    }
    else if(type == render_type::indexed)
    {
      numTriangles = state.num_triangles;
    }
    else
    {
      numTriangles = state.num_vertices-2;
    }
    printf("\t### num_triangles: %d\n",numTriangles);

    for(int triangleIndex = 0; triangleIndex < numTriangles; triangleIndex ++)
    {

      // printf("\t### traingleIndex: %d\n", triangleIndex);
      float x[state.num_vertices] = {};
      float y[state.num_vertices] = {};
      float z[state.num_vertices] = {};
      vec4 c[state.num_vertices] = {};
      float xMinObject = 1;
      float xMaxObject = -1;
      float yMinObject = 1;
      float yMaxObject = -1;

      for(int vertexNumber = 0; vertexNumber < 3; vertexNumber++)
      {
        int vertexOffset = 0;
        if(type == render_type::triangle)
        {
          vertexOffset = (triangleIndex*3 + vertexNumber)*state.floats_per_vertex;
        }
        else if(type == render_type::indexed)
        {
          vertexOffset = state.index_data[triangleIndex*3 + vertexNumber]*state.floats_per_vertex;
        }
        else if(type == render_type::fan)
        {
          if(vertexNumber > 0)
          {
            vertexOffset = (triangleIndex + vertexNumber)*state.floats_per_vertex;
          }
          else
          {
            vertexOffset = 0;
          }
        }
        else if(type == render_type::strip)
        {
          vertexOffset = (triangleIndex + vertexNumber)*state.floats_per_vertex;
        }

        tempVertex.data = &state.vertex_data[vertexOffset];
        state.vertex_shader(tempVertex, tempGeometry, state.uniform_data);
        x[vertexNumber] = tempGeometry.gl_Position.x[0]/tempGeometry.gl_Position.x[3];
        y[vertexNumber] = tempGeometry.gl_Position.x[1]/tempGeometry.gl_Position.x[3];
        z[vertexNumber] = tempGeometry.gl_Position.x[2]/tempGeometry.gl_Position.x[3];
        for(int fvi = 0; fvi < MAX_FLOATS_PER_VERTEX; fvi++)
        {
          if(state.interp_rules[fvi] == interp_type::flat)
          {
            tempFragment.data[fvi] = state.vertex_data[triangleIndex*state.floats_per_vertex + fvi];
          }
          else
          {
            tempFragment.data[fvi] = state.vertex_data[vertexOffset + fvi];
          }
        }
        state.fragment_shader(tempFragment, tempOutput, state.uniform_data);
        c[vertexNumber] = tempOutput.output_color;

        xMinObject = std::min(xMinObject, x[vertexNumber]);
        xMaxObject = std::max(xMaxObject, x[vertexNumber]);
        yMinObject = std::min(yMinObject, y[vertexNumber]);
        yMaxObject = std::max(yMaxObject, y[vertexNumber]);
      }

      state.num_triangles = state.num_vertices/3;








      // printf("colors0: %f, %f, %f, %f\n",c[0].x[0], c[0].x[1], c[0].x[2], c[0].x[3]);
      // printf("colors1: %f, %f, %f, %f\n",c[1].x[0], c[1].x[1], c[1].x[2], c[1].x[3]);
      // printf("colors2: %f, %f, %f, %f\n",c[2].x[0], c[2].x[1], c[2].x[2], c[2].x[3]);


      int xMinPixel, xMaxPixel, yMinPixel, yMaxPixel;
      xMinPixel = std::max((float)0,(float)round(((xMinObject+1)*state.image_width)/2));
      xMaxPixel = std::min((float)state.image_width,(float)round(((xMaxObject+1)*state.image_width)/2));
      yMinPixel = std::max((float)0,(float)round(((yMinObject+1)*state.image_height)/2));
      yMaxPixel = std::min((float)state.image_height,(float)round(((yMaxObject+1)*state.image_height)/2));

      vec3 A, B, C, P;
      float alpha, beta, gamma, areaABC;
      A = vec3(x[0], y[0], 0);
      B = vec3(x[1], y[1], 0);
      C = vec3(x[2], y[2], 0);
      // printf("A: %f,%f\tB: %f,%f\tC: %f,%f\n",A.x[0],A.x[1],B.x[0],B.x[1],C.x[0],C.x[1]);
      areaABC = cross(B - C, A - C).magnitude();
      // printf("areaABC: %f\n", areaABC);
      vec3 pixelSize = {2/(float)state.image_width, 2/(float)state.image_height,0};
      vec3 origin, e1, e2;
      float k0Alpha, k1Alpha, k2Alpha, k0Beta, k1Beta, k2Beta, k0Gamma, k1Gamma, k2Gamma;
      origin = vec3(-1,-1,0)+pixelSize/(float)2;
      e1 = origin + pixelSize*vec3(1,0,0);
      e2 = origin + pixelSize*vec3(0,1,0);

      k0Alpha = cross(C - B, origin - B).magnitude()/areaABC;// - xMinObject;
      k1Alpha = cross(C - B, e1 - B).magnitude()/areaABC - k0Alpha;
      k2Alpha = cross(C - B, e2 - B).magnitude()/areaABC - k0Alpha;
      // printf("k0Alpha: %f, k1Alpha: %f, k2Alpha: %f\n",k0Alpha, k1Alpha, k2Alpha);

      k0Beta = cross(A - C, origin - C).magnitude()/areaABC;// - xMinObject;
      k1Beta = cross(A - C, e1 - C).magnitude()/areaABC - k0Beta;
      k2Beta = cross(A - C, e2 - C).magnitude()/areaABC - k0Beta;
      // printf("k0Beta: %f, k1Beta: %f, k2Beta: %f\n", k0Beta, k1Beta, k2Beta);

      k0Gamma = cross(B - A, origin - A).magnitude()/areaABC;// - xMinObject;
      k1Gamma = cross(B - A, e1 - A).magnitude()/areaABC - k0Gamma;
      k2Gamma = cross(B - A, e2 - A).magnitude()/areaABC - k0Gamma;
      // printf("k0Gamma: %f, k1Gamma: %f, k2Gamma: %f\n", k0Gamma, k1Gamma, k2Gamma);

      P.x[0] = xMinPixel*pixelSize.x[0]-1+(pixelSize.x[0]/(float)2);
      P.x[1] = yMinPixel*pixelSize.x[1]-1+(pixelSize.x[1]/(float)2);
      P.x[2] = 0;
      alpha = cross(B - C, P - C).magnitude()/areaABC;
      beta = cross(A - C, P - C).magnitude()/areaABC;
      gamma = cross(A - B, P - B).magnitude()/areaABC;

      // alpha = k0Alpha + k1Alpha*xMinPixel + k2Alpha*yMinPixel;
      // beta =  k0Beta  + k1Beta *xMinPixel + k2Beta *yMinPixel;
      // gamma = k0Gamma + k1Gamma*xMinPixel + k2Gamma*yMinPixel);

      int imageIndex = yMinPixel*state.image_width + xMinPixel;
      int subWidthPixel = (xMaxPixel - xMinPixel);
      float subWidthObject = (xMaxObject - xMinObject);
      // float subWidthObject = subWidthPixel*k1Alpha;
      for(int indexY = yMinPixel; indexY <= yMaxPixel; indexY++)
      {
        for(int indexX = xMinPixel; indexX <= xMaxPixel; indexX++)
        {
      // imageIndex = 0;
      // for(int indexY = 0; indexY < state.image_height; indexY++)
      // {
      //   for(int indexX = 0; indexX < state.image_width; indexX++)
      //   {
          P.x[0] = indexX*pixelSize.x[0]-1+(pixelSize.x[0]/(float)2);
          P.x[1] = indexY*pixelSize.x[1]-1+(pixelSize.x[1]/(float)2);
          P.x[2] = 0;
          float alphaR = cross(C - B, P - B).magnitude()/areaABC;
          float betaR = cross(A - C, P - C).magnitude()/areaABC;
          float gammaR = cross(B - A, P - A).magnitude()/areaABC;
          // printf("xPixel: %d, yPixel: %d\n",indexX, indexY);
          // printf("Px,y: %f,%f\t alpha %f, beta %f, gamma %f, sum %f\n", P.x[0], P.x[1], alpha,beta,gamma, alpha+beta+gamma - 0.004167);
          // //
          // printf("REAL\t\t\t alpha %f, beta %f, gamma %f, sum %f\n", alphaR,betaR,gammaR, alphaR+betaR+gammaR);
          // exit(0);
          // printf("indexX: %d, indexY: %d, alpha: %f, beta: %f, gamma: %f, sum: %f\n",indexX,indexY, alpha, beta, gamma, alpha+beta+gamma);
          if(alphaR >= -1e-4 && betaR >= -1e-4 && gammaR >= -1e-4 && std::abs(1 - (alphaR + betaR + gammaR)) < 3e-4)
          {
            float bz = alphaR*z[0] + betaR*z[1] + gammaR*z[2];
            // printf("z old:
            if(bz < state.image_depth[imageIndex] && bz >= -1 && bz <= 1)
            {
              vec4 bc = alphaR*c[0] + betaR*c[1] + gammaR*c[2];
              // printf("colors: %f, %f, %f, %f\n",c[1].x[0], c[1].x[1], c[1].x[2], c[1].x[3]);
              bc *= 255;
              bc = componentwise_min(bc,vec4(255,255,255,0));
              state.image_color[imageIndex] = make_pixel((int)bc.x[0],(int)bc.x[1],(int)bc.x[2]);
              state.image_depth[imageIndex] = bz;
            }
          }
          alpha += k1Alpha;
          beta  += k1Beta;
          gamma += k1Gamma;
          imageIndex += 1;

        }
        // if(indexY == yMinPixel+1)
        // {
        //   exit(0);
        // }
        alpha += k2Alpha - (subWidthPixel+1)*k1Alpha;
        beta += k2Beta - (subWidthPixel+1)*k1Beta;
        gamma += k2Gamma - (subWidthPixel+1)*k1Gamma;
        imageIndex += state.image_width - subWidthPixel - 1;
      }

    }
    delete [] tempGeometry.data;
    delete [] tempFragment.data;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    std::cout<<"TODO: implement rasterization"<<std::endl;
}
