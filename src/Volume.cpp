#include <iostream>
#include <fstream>
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
    if (vertex_shader_file.empty() || fragment_shader_file.empty())
    {
        return;
    }
    // create vertex shader
    vertexShaderID = glCreateShader(GL_VERTEX_SHADER);
    // create fragment shader
    fragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);
    // read src
    std::string vertex_source = ReadFile(vertex_shader_file);
    std::string fragment_source = ReadFile(fragment_shader_file);
    const char *v_source = vertex_source.c_str();
    const char *f_source = fragment_source.c_str();
    if (vertex_source.empty() || fragment_source.empty())
    {
        return;
    }
    // bind shader
    glShaderSource(vertexShaderID, 1, &v_source, 0);
    glShaderSource(fragmentShaderID, 1, &f_source, 0);

    glCompileShader(vertexShaderID);
    glCompileShader(fragmentShaderID);

    programID = glCreateProgram();

    glAttachShader(programID, vertexShaderID);
    glAttachShader(programID, fragmentShaderID);

    glLinkProgram(programID);

    /* get a handle for uniform */
    volumeTexID = createPyroclasticVolume(1, radius);
    cameraPosID = glGetUniformLocation(programID, "eyePos");
    LightPosID = glGetUniformLocation(programID, "lightPos");
    LightIntensityID = glGetUniformLocation(programID, "lightIntensity");
    absorptionID = glGetUniformLocation(programID, "absorption");
    MatrixID = glGetUniformLocation(programID, "MVP");

    lightPos = glm::vec3(0, 0, 4);
    lightIntensity = glm::vec3(1.0f, 1.0f, 1.0f);
    absorption = 0.1f;
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

    glDetachShader(programID, vertexShaderID);
    glDetachShader(programID, fragmentShaderID);

    glDeleteShader(vertexShaderID);
    glDeleteShader(fragmentShaderID);

    glDeleteProgram(programID);
}

void Volume::update()
{
    m_camera.update();
    cameraPos = m_camera.getPos();
    lightPos = glm::vec3(0, 0, 4);
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
    glDisableVertexAttribArray(1);
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

GLuint Volume::createPyroclasticVolume(int n, float r)
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

    float frequency = 3.0f / n;
    float center = n / 2.0f + 0.5f;

    for (int x = 0; x < n; x++)
    {
        for (int y = 0; y < n; ++y)
        {
            for (int z = 0; z < n; ++z)
            {
                float dx = center - x;
                float dy = center - y;
                float dz = center - z;

                float off = fabsf(glm::perlin(glm::vec3(x, y, z), glm::vec3(frequency)));

                float d = sqrtf(dx * dx + dy * dy + dz * dz) / (n);

                *ptr++ = ((d - off) < r) ? 255 : 0;
            }
        }
    }

    // upload
    glTexImage3D(target,
                 0,
                 GL_LUMINANCE,
                 n,
                 n,
                 n,
                 0,
                 GL_LUMINANCE,
                 GL_UNSIGNED_BYTE,
                 data);

    glBindTexture(target, 0);
    delete[] data;

    return texid;
}