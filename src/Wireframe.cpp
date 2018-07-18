#include <fstream>
#include <sstream>
#include <vector>
#include "Wireframe.hpp"
#include "constants.hpp"

Wireframe::Wireframe()
{
    initialize();
}

Wireframe::~Wireframe()
{
    glDeleteBuffers(1, &indexID);
    glDeleteBuffers(1, &vboID);
    glDeleteVertexArrays(1, &vaoID);

    glDeleteProgram(programID);
}

void Wireframe::update()
{
}

void Wireframe::draw() const
{
    // Draw the triangles !
    glDrawElements(
        GL_LINES,        // mode
        24,              // count
        GL_UNSIGNED_INT, // type
        (void *)0        // element array buffer offset
    );
}

GLuint Wireframe::getProgramID() const
{
    return programID;
}
GLuint Wireframe::getVaoID() const
{
    return vaoID;
}
GLuint Wireframe::getMatrixID() const
{
    return MatrixID;
}

/* private */

void Wireframe::initialize()
{
    initVAO();
    initShaders();
}

void Wireframe::initVAO()
{
    /* domain cube */
    float vertices[8][3] = {
        {-1.0f, -1.0f, -1.0f},
        {1.0f, -1.0f, -1.0f},
        {-1.0f, 1.0f, -1.0f},
        {-1.0f, -1.0f, 1.0f},
        {1.0f, 1.0f, -1.0f},
        {-1.0f, 1.0f, 1.0f},
        {1.0f, -1.0f, 1.0f},
        {1.0f, 1.0f, 1.0f}};

    const unsigned int indices[12][2] = {
        {0, 1}, {1, 6}, {6, 3}, {3, 0}, {0, 2}, {1, 4}, {6, 7}, {3, 5}, {2, 4}, {4, 7}, {7, 5}, {5, 2}};

    for (int i = 0; i < 8; ++i)
    {
        vertices[i][0] *= (float)ratio[0];
        vertices[i][1] *= (float)ratio[1];
        vertices[i][2] *= (float)ratio[2];
    }

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

void Wireframe::initShaders()
{
    std::string vertex_shader_file = std::string("./src/shader/wireframe.vert");
    std::string fragment_shader_file = std::string("./src/shader/wireframe.frag");

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

    MatrixID = glGetUniformLocation(programID, "MVP");
}

std::string Wireframe::ReadFile(const std::string &filename)
{
    std::ifstream ifs(filename);
    if (ifs.is_open())
        return nullptr;
    std::istreambuf_iterator<char> ifs_begin(ifs);
    std::istreambuf_iterator<char> ifs_end;
    std::string file_string(ifs_begin, ifs_end);
    return file_string;
}
