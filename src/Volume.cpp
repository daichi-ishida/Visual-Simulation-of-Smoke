#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "Volume.hpp"
#include "constants.hpp"

Volume::Volume(MACGrid *grids, const std::string &vertex_shader_file, const std::string &fragment_shader_file) : m_grids(grids), m_camera()
{
    /* domain cube */
    const float vertices[8][3] = {
        {-1.0f, -1.0f, -1.0f},
        {1.0f, -1.0f, -1.0f},
        {-1.0f, 1.0f, -1.0f},
        {-1.0f, -1.0f, 1.0f},
        {1.0f, 1.0f, -1.0f},
        {-1.0f, 1.0f, 1.0f},
        {1.0f, -1.0f, 1.0f},
        {1.0f, 1.0f, 1.0f}};

    const unsigned int indices[12][3] = {
        {3, 5, 7}, {3, 7, 6}, {1, 6, 7}, {1, 7, 4}, {0, 1, 4}, {0, 4, 2}, {0, 2, 5}, {0, 5, 3}, {2, 5, 7}, {2, 7, 4}, {0, 1, 6}, {0, 6, 3}};

    // generate VAO
    glGenVertexArrays(1, &m_vertex_array_object);
    // set current VAO
    glBindVertexArray(m_vertex_array_object);

    // generate VBO
    glGenBuffers(1, &m_vertex_buffer_object);
    glBindBuffer(GL_ARRAY_BUFFER, m_vertex_buffer_object);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    // generate index buffer
    glGenBuffers(1, &m_index_buffer_object);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_index_buffer_object);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

    glBindVertexArray(0);

    /* shader */
    // create shaders
    vertexShaderID = glCreateShader(GL_VERTEX_SHADER);
    fragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);

    // Read the Vertex Shader code from the file
    std::string vertexShaderCode;
    std::ifstream vertexShaderStream(vertex_shader_file, std::ios::in);
    if (vertexShaderStream.is_open())
    {
        std::stringstream sstr;
        sstr << vertexShaderStream.rdbuf();
        vertexShaderCode = sstr.str();
        vertexShaderStream.close();
    }
    else
    {
        printf("Impossible to open %s. Are you in the right directory ? Don't forget to read the FAQ !\n", vertex_shader_file);
        getchar();
        return;
    }

    // Read the Fragment Shader code from the file
    std::string fragmentShaderCode;
    std::ifstream fragmentShaderStream(fragment_shader_file, std::ios::in);
    if (fragmentShaderStream.is_open())
    {
        std::stringstream sstr;
        sstr << fragmentShaderStream.rdbuf();
        fragmentShaderCode = sstr.str();
        fragmentShaderStream.close();
    }

    GLint Result = GL_FALSE;
    int InfoLogLength;

    // Compile Vertex Shader
    printf("Compiling shader : %s\n", vertex_shader_file);
    char const *vertexSourcePointer = vertexShaderCode.c_str();
    glShaderSource(vertexShaderID, 1, &vertexSourcePointer, NULL);
    glCompileShader(vertexShaderID);

    // Check Vertex Shader
    glGetShaderiv(vertexShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(vertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    if (InfoLogLength > 0)
    {
        std::vector<char> vertexShaderErrorMessage(InfoLogLength + 1);
        glGetShaderInfoLog(vertexShaderID, InfoLogLength, NULL, &vertexShaderErrorMessage[0]);
        printf("%s\n", &vertexShaderErrorMessage[0]);
    }

    // Compile Fragment Shader
    printf("Compiling shader : %s\n", fragment_shader_file);
    char const *fragmentSourcePointer = fragmentShaderCode.c_str();
    glShaderSource(fragmentShaderID, 1, &fragmentSourcePointer, NULL);
    glCompileShader(fragmentShaderID);

    // Check Fragment Shader
    glGetShaderiv(fragmentShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(fragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    if (InfoLogLength > 0)
    {
        std::vector<char> fragmentShaderErrorMessage(InfoLogLength + 1);
        glGetShaderInfoLog(fragmentShaderID, InfoLogLength, NULL, &fragmentShaderErrorMessage[0]);
        printf("%s\n", &fragmentShaderErrorMessage[0]);
    }

    // Link the program
    printf("Linking program\n");
    programID = glCreateProgram();
    glAttachShader(programID, vertexShaderID);
    glAttachShader(programID, fragmentShaderID);
    glLinkProgram(programID);

    // Check the program
    glGetProgramiv(programID, GL_LINK_STATUS, &Result);
    glGetProgramiv(programID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    if (InfoLogLength > 0)
    {
        std::vector<char> programErrorMessage(InfoLogLength + 1);
        glGetProgramInfoLog(programID, InfoLogLength, NULL, &programErrorMessage[0]);
        printf("%s\n", &programErrorMessage[0]);
    }

    glDetachShader(programID, vertexShaderID);
    glDetachShader(programID, fragmentShaderID);

    glDeleteShader(vertexShaderID);
    glDeleteShader(fragmentShaderID);

    /* get a handle for uniform */
    volumeTexID = createPyroclasticVolume(Nx, RADIUS_THRESHOLD);
    cameraPosID = glGetUniformLocation(programID, "eyePos");
    LightPosID = glGetUniformLocation(programID, "lightPos");
    LightIntensityID = glGetUniformLocation(programID, "lightIntensity");
    absorptionID = glGetUniformLocation(programID, "absorption");
    MatrixID = glGetUniformLocation(programID, "MVP");

    lightPos = glm::vec3(0, 4, 4);
    lightIntensity = glm::vec3(1.0f, 1.0f, 1.0f);
    absorption = 0.8f;
    glm::mat4 ProjectionMatrix = m_camera.getProjectionMat();
    glm::mat4 ViewMatrix = m_camera.getViewMat();
    glm::mat4 ModelMatrix = glm::mat4(1.0);
    MVP = ProjectionMatrix * ViewMatrix * ModelMatrix;
}

Volume::~Volume()
{
    glDeleteBuffers(1, &m_index_buffer_object);
    glDeleteBuffers(1, &m_vertex_buffer_object);
    glDeleteVertexArrays(1, &m_vertex_array_object);

    glDeleteProgram(programID);
}

void Volume::update()
{
    m_camera.update();
    cameraPos = m_camera.getPos();
    lightPos = glm::vec3(0, 4, 4);
    lightIntensity = glm::vec3(1.0f, 1.0f, 1.0f);

    // matrix
    glm::mat4 ProjectionMatrix = m_camera.getProjectionMat();
    glm::mat4 ViewMatrix = m_camera.getViewMat();
    glm::mat4 ModelMatrix = glm::mat4(1.0);
    MVP = ProjectionMatrix * ViewMatrix * ModelMatrix;
}

void Volume::draw()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Use our shader
    glUseProgram(programID);

    // set uniform
    glUniform3f(cameraPosID, cameraPos.x, cameraPos.y, cameraPos.z);
    glUniform3f(LightPosID, lightPos.x, lightPos.y, lightPos.z);
    glUniform3f(LightIntensityID, lightIntensity.x, lightIntensity.y, lightIntensity.z);
    glUniform1f(absorptionID, absorption);
    glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);

    // Bind volume texture in Texture Unit 0
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, volumeTexID);
    // Set "VolumeTextureSampler" sampler to user Texture Unit 0
    glUniform1i(volumeTexID, 0);

    // 1rst attribute buffer : vertices
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, m_vertex_buffer_object);
    glVertexAttribPointer(
        0,        // attribute
        3,        // size
        GL_FLOAT, // type
        GL_FALSE, // normalized?
        0,        // stride
        (void *)0 // array buffer offset
    );
    // // 2nd attribute buffer : texcoords
    // glEnableVertexAttribArray(0);
    // glBindBuffer(GL_ARRAY_BUFFER, m_texture_buffer_object);
    // glVertexAttribPointer(
    //     1,        // attribute
    //     3,        // size
    //     GL_FLOAT, // type
    //     GL_FALSE, // normalized?
    //     0,        // stride
    //     (void *)0 // array buffer offset
    // );

    // Index buffer
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_index_buffer_object);

    // Draw the triangles !
    glDrawElements(
        GL_TRIANGLES,    // mode
        36,              // count
        GL_UNSIGNED_INT, // type
        (void *)0        // element array buffer offset
    );

    glDisableVertexAttribArray(0);
    glBindTexture(GL_TEXTURE_3D, 0);
    glUseProgram(0);
}

/* private */

std::string Volume::ReadFile(const std::string &filename)
{
    std::ifstream ifs(filename);
    if (ifs.is_open())
        return nullptr;
    std::istreambuf_iterator<char> ifs_begin(ifs);
    std::istreambuf_iterator<char> ifs_end;
    std::string file_string(ifs_begin, ifs_end);
    return file_string;
}

GLuint Volume::createPyroclasticVolume(int n, float r_th)
{

    GLuint texid;
    glGenTextures(1, &texid);

    GLenum target = GL_TEXTURE_3D;
    GLenum filter = GL_LINEAR;
    GLenum address = GL_CLAMP_TO_BORDER;

    glBindTexture(target, texid);

    glTexParameteri(target, GL_TEXTURE_MAG_FILTER, filter);
    glTexParameteri(target, GL_TEXTURE_MIN_FILTER, filter);

    glTexParameteri(target, GL_TEXTURE_WRAP_S, address);
    glTexParameteri(target, GL_TEXTURE_WRAP_T, address);
    glTexParameteri(target, GL_TEXTURE_WRAP_R, address);

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    GLbyte *data = new GLbyte[SIZE];
    GLbyte *ptr = data;

    glm::vec3 frequency = glm::vec3(3.0f / Nx, 3.0f / Ny, 3.0f / Nz);
    glm::vec3 center = 0.5f * glm::vec3(Nx, Ny, Nz);

    for (int x = 0; x < Nx; ++x)
    {
        for (int y = 0; y < Ny; ++y)
        {
            for (int z = 0; z < Nz; ++z)
            {
                glm::vec3 r = glm::vec3(x, y, z);
                glm::vec3 dr = center - r;

                float off = fabsf(glm::perlin(r, frequency));
                // float off = 0.0;

                float d = glm::length(dr) / Nx;

                *ptr++ = ((d - off) < r_th) ? 255 : 0;
            }
        }
    }

    // upload
    glTexImage3D(target,
                 0,
                 GL_LUMINANCE,
                 Nx,
                 Ny,
                 Nz,
                 0,
                 GL_LUMINANCE,
                 GL_UNSIGNED_BYTE,
                 data);

    glBindTexture(target, 0);
    delete[] data;

    return texid;
}