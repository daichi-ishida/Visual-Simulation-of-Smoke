#include <fstream>
#include <sstream>
#include <vector>
#include "Volume.hpp"
#include "constants.hpp"

Volume::Volume(std::shared_ptr<MACGrid> grids) : m_grids(grids)
{
    initialize();
}

Volume::~Volume()
{
    glDeleteBuffers(1, &indexID);
    glDeleteBuffers(1, &vboID);
    glDeleteVertexArrays(1, &vaoID);

    glDeleteProgram(programID);
}

void Volume::update()
{
    GLenum target = GL_TEXTURE_3D;
    GLenum filter = GL_LINEAR;
    GLenum address = GL_CLAMP_TO_BORDER;
    glBindTexture(target, volumeTexID);

    glTexParameteri(target, GL_TEXTURE_MAG_FILTER, filter);
    glTexParameteri(target, GL_TEXTURE_MIN_FILTER, filter);

    glTexParameteri(target, GL_TEXTURE_WRAP_S, address);
    glTexParameteri(target, GL_TEXTURE_WRAP_T, address);
    glTexParameteri(target, GL_TEXTURE_WRAP_R, address);

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    GLubyte *data = new GLubyte[SIZE];
    GLubyte *ptr = data;

    for (int z = 0; z < Nz; ++z)
    {
        for (int y = 0; y < Ny; ++y)
        {
            for (int x = 0; x < Nx; ++x)
            {
                float f = (float)m_grids->density(x, y, z);
                *ptr++ = std::max(0, std::min(255, (int)std::floor(f * 256.0)));
            }
        }
    }
    glTexImage3D(target,
                 0,
                 GL_RED,
                 Nx,
                 Ny,
                 Nz,
                 0,
                 GL_RED,
                 GL_UNSIGNED_BYTE,
                 data);

    glBindTexture(target, 0);
}

void Volume::draw() const
{
    glUniform1f(absorptionID, ABSORPTION);
    glUniform1f(numID, 2 * N);
    glUniform3f(ratioID, (float)ratio[0], (float)ratio[1], (float)ratio[2]);

    // Bind volume texture in Texture Unit 0
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, volumeTexID);

    // Draw the triangles !
    glDrawElements(
        GL_TRIANGLES,    // mode
        36,              // count
        GL_UNSIGNED_INT, // type
        (void *)0        // element array buffer offset
    );

    glBindTexture(GL_TEXTURE_3D, 0);
}

/* getID for Uniform */

GLuint Volume::getProgramID() const
{
    return programID;
}
GLuint Volume::getVaoID() const
{
    return vaoID;
}
GLuint Volume::getCamPosID() const
{
    return cameraPosID;
}
GLuint Volume::getLightPosID() const
{
    return LightPosID;
}
GLuint Volume::getLightIntensityID() const
{
    return LightIntensityID;
}
GLuint Volume::getMatrixID() const
{
    return MatrixID;
}

/* private */

void Volume::initialize()
{
    initVAO();
    initShaders();
}

void Volume::initVAO()
{
    /* domain cube */
    const float vertices[8][3] = {
        {0.0f, 0.0f, 0.0f},
        {1.0f, 0.0f, 0.0f},
        {0.0f, 1.0f, 0.0f},
        {0.0f, 0.0f, 1.0f},
        {1.0f, 1.0f, 0.0f},
        {0.0f, 1.0f, 1.0f},
        {1.0f, 0.0f, 1.0f},
        {1.0f, 1.0f, 1.0f}};

    const unsigned int indices[12][3] = {
        {3, 5, 7}, {3, 7, 6}, {1, 6, 7}, {1, 7, 4}, {0, 1, 4}, {0, 4, 2}, {0, 2, 5}, {0, 5, 3}, {2, 5, 7}, {2, 7, 4}, {0, 1, 6}, {0, 6, 3}};

    // generate VAO
    glGenVertexArrays(1, &vaoID);
    // set current VAO
    glBindVertexArray(vaoID);

    // generate VBO
    glGenBuffers(1, &vboID);
    glBindBuffer(GL_ARRAY_BUFFER, vboID);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(
        0,        // attribute
        3,        // size
        GL_FLOAT, // type
        GL_FALSE, // normalized?
        0,        // stride
        (void *)0 // array buffer offset
    );

    // generate index buffer
    glGenBuffers(1, &indexID);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexID);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

    glBindVertexArray(0);
}

void Volume::initShaders()
{
    std::string vertex_shader_file = std::string("./src/shader/volume.vert");
    std::string fragment_shader_file = std::string("./src/shader/volume.frag");

    GLuint vertexShaderID;
    GLuint fragmentShaderID;
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
        printf("Impossible to open %s. Are you in the right directory ? Don't forget to read the FAQ !\n", vertex_shader_file.c_str());
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
    else
    {
        printf("Impossible to open %s. Are you in the right directory ? Don't forget to read the FAQ !\n", fragment_shader_file.c_str());
        getchar();
        return;
    }

    GLint Result = GL_FALSE;
    int InfoLogLength;

    // Compile Vertex Shader
    printf("Compiling shader : %s\n", vertex_shader_file.c_str());
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
    printf("Compiling shader : %s\n", fragment_shader_file.c_str());
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
    cameraPosID = glGetUniformLocation(programID, "eyePos");
    LightPosID = glGetUniformLocation(programID, "lightPos");
    LightIntensityID = glGetUniformLocation(programID, "lightIntensity");
    MatrixID = glGetUniformLocation(programID, "MVP");
    absorptionID = glGetUniformLocation(programID, "absorption");
    numID = glGetUniformLocation(programID, "num");
    ratioID = glGetUniformLocation(programID, "ratio");
}

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