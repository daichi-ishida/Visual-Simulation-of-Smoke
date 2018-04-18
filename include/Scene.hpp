#pragma once

class Scene
{
  public:
    Scene();
    ~Scene();

    void writeData();

  private:
    void writeData_inVtuFormat();

    int m_file_num;
};
