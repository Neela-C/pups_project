
#include <glad/glad.h>

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <limits>
#include <functional>

#include "Geometry.h"
#include "GLDebug.h"
#include "Log.h"
#include "ShaderProgram.h"
#include "Shader.h"
#include "Texture.h"
#include "Window.h"

#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"


#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
//#include "Scene.h"
//#include "Lighting.h"
using namespace std;
int view_index = 0;
// color map
glm::vec3 color1 = glm::vec3(0.0, 0.0, 0.5);   // dark blue
glm::vec3 color2 = glm::vec3(0.0, 0.5, 1.0);   // light blue
glm::vec3 color5 = glm::vec3(0.9, 0.9, 0.0);   // yellow
glm::vec3 color3 = glm::vec3(0.5, 0.0, 0.0);   // dark red
glm::vec3 color4 = glm::vec3(1.0, 0.0, 0.0);   // bright red

glm::vec3 lightPosition(0.0f, 0.0f, 0.0f);
//glm::vec3 Pointlight_position = glm::vec3(0.0f, 0.0f, -10.0f);


float threshold1 = 0.0;
float threshold2 = 0.25;
float threshold3 = 0.5;
float threshold4 = 0.75;
float threshold5 = 1.0;


/*
std::vector<glm::vec3> get_normals(std::vector<glm::vec3> points) {
	std::vector<glm::vec3> normals(points.size(), glm::vec3(0.0f));

	for (size_t i = 0; i < points.size()-3; i += 3) {
		glm::vec3 edge1 = points[i + 1] - points[i];
		glm::vec3 edge2 = points[i + 2] - points[i];
		glm::vec3 normal = glm::normalize(glm::cross(edge1, edge2));
		normals[i] += normal;
		normals[i + 1] += normal;
		normals[i + 2] += normal;
	}

	for (size_t i = 0; i < normals.size(); i++) {
		normals[i] = glm::normalize(normals[i]);
	}
	return normals;
}*/

std::vector<glm::vec3> get_normals(std::vector<glm::vec3> points) {
	std::vector<glm::vec3> normals(points.size(), glm::vec3(0.0f));
	std::vector<int> counts(points.size(), 0);

	for (size_t i = 0; i < points.size(); i += 3) {
		glm::vec3 v1 = points[i];
		glm::vec3 v2 = points[i + 1];
		glm::vec3 v3 = points[i + 2];
		glm::vec3 edge1 = v2 - v1;
		glm::vec3 edge2 = v3 - v1;
		glm::vec3 normal = glm::normalize(glm::cross(edge1, edge2));
		normals[i] += normal;
		normals[i + 1] += normal;
		normals[i + 2] += normal;
		counts[i]++;
		counts[i + 1]++;
		counts[i + 2]++;
	}

	for (size_t i = 0; i < normals.size(); i++) {
		if (counts[i] > 0) {
			normals[i] /= counts[i];
		}
		if (glm::length(normals[i]) > 0) {
			normals[i] = glm::normalize(normals[i]);
		}
	}

	return normals;
}



void print_vec3(glm::vec3 vec) {
	std::cerr << "{ ";
	std::cerr << vec[0] << ", " << vec[1] << ", " << vec[2];
	std::cerr << " }, ";
}

void print_vec(std::vector<glm::vec3> vec) {
	for (auto i : vec) {
		std::cerr << "{ ";
		print_vec3(i);
		std::cerr << " }\n";
	}
}

void print_mat4(glm::mat4 vec) {
	std::cerr << "[ ";
	for (auto i = 0; i < 4; i++) {
		for (auto j = 0; j < 4; j++) {
			std::cout << vec[i][j] << ", ";
		}
		std::cout << std::endl;
	}
	std::cerr << " ]\n";
}


/*
 * https://learnopengl.com/Getting-started/Camera - referenced
 */
float width = 800, height = 800;
glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f, 2.4f);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);
const float cameraSpeed = 0.3f;

bool firstMouse = true;
float yaw = -90.0f;
float pitch = 0.0f;
float lastX = width / 2.0;
float lastY = height / 2.0;

struct State {
    bool rightMousePressed = false;
    bool leftMousePressed = false;
    bool leftMouseReleased = false;
    bool leftMouseHeld = false;
    bool cameraForward = false;
    bool cameraBackward = false;
    bool cameraLeft = false;
    bool cameraRight = false;
    bool cameraUp = false;
    bool cameraDown = false;
	bool cameraReset = false;

	bool controlUp = false;
	bool controlDown = false;

    glm::vec2 mouseCoord;
    glm::vec2 mouseDownCoord;
    glm::vec2 mouseUpCoord;

    bool mode = 1;
    bool revolve = 0;
    bool wire = 0;
	bool pups = 1;
    int scene = 7;
};

void updateGPUGeometry(GPU_Geometry &gpuGeom, CPU_Geometry const &cpuGeom) {
    gpuGeom.bind();
    gpuGeom.setVerts(cpuGeom.verts);
    gpuGeom.setCols(cpuGeom.cols);
	gpuGeom.setNormals(cpuGeom.normals);
}

void updateHeightGPUGeometry(GPU_Geometry_height& gpuGeom, CPU_Geometry_height const& cpuGeom) {
	gpuGeom.bind();
	gpuGeom.setVerts(cpuGeom.verts);
	gpuGeom.setHeight(cpuGeom.heights);
	gpuGeom.setCols(cpuGeom.cols);
}

class MyCallbacks : public CallbackInterface {

public:
    MyCallbacks(int screenWidth, int screenHeight) :
            screenDimensions(screenWidth, screenHeight) {
    }

    virtual void keyCallback(int key, int scancode, int action, int mods) {
        if (key == GLFW_KEY_R && action == GLFW_PRESS) {
            state.scene = -1;
        }
        /*if (key == GLFW_KEY_1 && action == GLFW_PRESS) {
            state.scene = 1;
        }*/
        if (key == GLFW_KEY_2 && action == GLFW_PRESS) {
            state.scene = 2;
        }
        /*if (key == GLFW_KEY_3 && state.mode && action == GLFW_PRESS) {
            if (state.revolve == 0) {
                state.revolve = 1;
            } else {
                state.revolve = 0;
            }
        }*/
        /*if (key == GLFW_KEY_4 && action == GLFW_PRESS) {
            state.scene = 4;
        }
        if (key == GLFW_KEY_5 && action == GLFW_PRESS) {
            state.scene = 5;
        }
        if (key == GLFW_KEY_6 && action == GLFW_PRESS) {
            state.scene = 6;
        }*/
        if (key == GLFW_KEY_7 && action == GLFW_PRESS) {
            state.scene = 7;
			//state.mode = 0;
			state.pups = 1;
        }
        if (key == GLFW_KEY_TAB && action == GLFW_PRESS) {
			
            /*if (state.mode == 0) {
                state.mode = 1;
            } else {
                state.mode = 0;
            }*/
			state.mode = 1;
        }
        if (key == GLFW_KEY_CAPS_LOCK && state.mode && action == GLFW_PRESS) {
            if (state.wire == 0) {
                state.wire = 1;
            } else {
                state.wire = 0;
            }
        }
		if (key == GLFW_KEY_U && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
			state.controlUp = true;
		}
		if (key == GLFW_KEY_J && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
			state.controlDown = true;
		}
        if (key == GLFW_KEY_W && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
            state.cameraForward = true;
        }
        if (key == GLFW_KEY_S && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
            state.cameraBackward = true;
        }
        if (key == GLFW_KEY_A && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
            state.cameraLeft = true;
        }
        if (key == GLFW_KEY_D && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
            state.cameraRight = true;
        }
		if (key == GLFW_KEY_T && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
			state.cameraReset = true;
		}
        if (key == GLFW_KEY_UP && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
            state.cameraUp = true;
        }
        if (key == GLFW_KEY_DOWN && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
            state.cameraDown = true;
        }
        stateChanged = true;
    }

    virtual void mouseButtonCallback(int button, int action, int mods) {
        if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
            state.mouseDownCoord = state.mouseCoord;
            state.leftMouseHeld = true;
            state.leftMousePressed = true;
        } else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
            state.mouseUpCoord = state.mouseCoord;
            state.leftMouseHeld = false;
            state.leftMouseReleased = true;
        } else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
            state.mouseDownCoord = state.mouseCoord;
            state.rightMousePressed = true;
        }
        stateChanged = true;
    }

    virtual void cursorPosCallback(double xpos, double ypos) {
        xScreenPos = xpos;
        yScreenPos = ypos;
        state.mouseCoord = coordsGL();

        stateChanged = true;
    }

    virtual void scrollCallback(double xoffset, double yoffset) {
    }

    virtual void windowSizeCallback(int width, int height) {
        CallbackInterface::windowSizeCallback(width, height);
        screenDimensions.x = width;
        screenDimensions.y = height;
    }

    glm::vec2 coordsGL() {
        glm::vec2 startingVec(xScreenPos, yScreenPos);
        glm::vec2 shiftedVec = startingVec + glm::vec2(0.5f);
        glm::vec2 scaledToZeroOne = shiftedVec / glm::vec2(screenDimensions);

        glm::vec2 flippedY = glm::vec2(scaledToZeroOne.x, 1.0f - scaledToZeroOne.y);

        glm::vec2 final = flippedY * 2.0f - glm::vec2(1.0f);

        return final;
    }

    State &getState() {
        return state;
    }

    bool checkStateChanged() {
        return stateChanged;
    }

    void stateHandled() {
        stateChanged = false;
        state.rightMousePressed = false;
        state.leftMousePressed = false;
        state.leftMouseReleased = false;
        state.cameraForward = false;
        state.cameraBackward = false;
        state.cameraRight = false;
        state.cameraLeft = false;
        state.cameraUp = false;
        state.cameraDown = false;
		state.cameraReset = false;
    }

    void resetScene() {
        state.scene = 7;
        state.mode = 1;
        state.revolve = 0;
        state.wire = 0;
		state.pups = 1;
        cameraPos = glm::vec3(0.0f, 0.0f, 2.4f);
        cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
        cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);
		/*cameraPos = glm::vec3(0.0f, 0.75f, 0.75f);
		cameraFront = glm::vec3(0.0f, 0.5f, -0.5f);
		cameraUp = glm::vec3(0.0f, 0.5f, 0.5f);*/

    }

    glm::vec2 getMousePos() {
        return glm::vec2(xScreenPos, yScreenPos);
    }

private:
    State state;
    bool stateChanged = true;
    glm::ivec2 screenDimensions;
    double xScreenPos;
    double yScreenPos;
};


vector<glm::vec2> findCart(int arr1[], int arr2[], int n, int n1) ///cartesian product
{
	vector<glm::vec2> cart_product;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n1; j++)
			cart_product.push_back(glm::vec2(arr1[i], arr2[j]));
			//printf("{%d, %d}, ", arr1[i], arr2[j]);
	return cart_product;
}



glm::vec3 pups(vector<glm::vec3>& controlPoints, const float u) {
	
	glm::vec3 return_point(0.f,0.f,0.f);
	vector <float> W_i;
	float k = 10.f;
	float weight;
	int N = controlPoints.size(); ///number of control points
	float C = 3 / (float)N; ///////the number of knot values overlapping
	float sumOfWeights = 0;
	
	for (float i = 0; i < controlPoints.size(); i++) {
		//weight = exp(-k * pow((u - i / N), 2) / (C * C - pow((u - i / N), 2)));

		
		if (u - float(i) / N < C && u - float(i) / N > -C) { ////if its inside the interval
			weight = exp(-k * pow(u - i / N, 2) / (C * C - pow((u - i / N), 2)));

		}

		else {
			weight = 0;
		}
		
		//weight = exp(-k * pow((u - i / N), 2) / (C * C - pow((u - i / N), 2)));
		W_i.push_back(weight);
		
		sumOfWeights += weight;
		
	}
	
	for (int i = 0; i < controlPoints.size(); i++)
	{
		return_point += W_i[i] * controlPoints[i];
	//std:cout << "Im in!";
	}
	return_point = return_point / sumOfWeights;
	//E_p.push_back()
	return return_point;

}



float weight_f(int N, float u, float v, float i, float j, float C) {

	//float C = 10.f / (float)N; // 16
	float k = 2.f; // 2
	float weight;
	float weight_1;
	float weight_2;
	//weight = exp(-k * pow(d - i / N, 2) / (C * C - pow((d - i / N), 2)));
	u = u - i / N;
	v = v - j / N;
	//float C = 3.0f / (float)N;

	if (u <= C && u >= -C) { ////if its inside the interval

		if (v <= C && v >= -C) {


			weight_1 = exp(-k * pow(u, 2) / ((C * C) - (u * u)));
			weight_2 = exp(-k * pow(v, 2) / ((C * C) - (v * v)));
			weight = weight_1 * weight_2;
			
			
		}
		else {
			weight = 0;
		}
	}
	else {
		weight = 0;
	}

	
	return weight;
}


glm::vec3 pups_surface(vector<vector<glm::vec3>>& controlPoints, const float u, const float v) {
	vector <float> W_i;
	
	float weight;
	int N = controlPoints.size(); 	
	float C = 3.0f / (float)N; 	
	float sumOfWeights = 0;
	glm::vec3 weightedSumOfPoints(0, 0, 0);
	
	for (float i = 0; i < controlPoints.size(); i++) {
		for (float j = 0; j < controlPoints[i].size(); j++) {
			weight = weight_f(N, u, v, i, j, C);
			W_i.push_back(weight);
			weightedSumOfPoints += weight * controlPoints[i][j];
			sumOfWeights += weight;
		}
	}
	return weightedSumOfPoints / sumOfWeights;
}


void bezierCurve(vector<glm::vec3> &controlPoints, vector<glm::vec3> &bezierPoints, int degree) {
    for (float u = 0; u < 1; u += 0.01) {
        for (int i = 1; i < degree; i++) {
            for (int j = 0; j < (degree - i); j++) {
                controlPoints[j] = (1 - u) * controlPoints[j] + u * controlPoints[j + 1];
            }
        }
        if (controlPoints.size() > 0) {
            bezierPoints.push_back(controlPoints[0]);
        }
    }
}



void splineCurve(vector<glm::vec3> &controlPoints, vector<glm::vec3> &splinePoints, int iteration) {
    if (iteration > 0) {
        controlPoints.assign(splinePoints.begin(), splinePoints.end());
    }
    if (controlPoints.size() > 0) {
        splinePoints.resize((controlPoints.size() - 1) * 2);
    }
    int points = controlPoints.size();
    int j = 0;
    for (int i = 0; i < points - 1; i++) {
        if (j == 0) {
            splinePoints[j] = controlPoints[i];
            splinePoints[j + 1] = (0.5f) * controlPoints[i] + (0.5f) * controlPoints[i + 1];
        } else if (j == splinePoints.size() - 2) {
            splinePoints[j] = (0.5f) * controlPoints[i] + (0.5f) * controlPoints[i + 1];
            splinePoints[j + 1] = controlPoints[i + 1];
        } else {
            splinePoints[j] = (0.75f) * controlPoints[i] + (0.25f) * controlPoints[i + 1];
            splinePoints[j + 1] = (0.25f) * controlPoints[i] + (0.75f) * controlPoints[i + 1];
			
        }

        j += 2;
    }
}


void reorderSorCurve(vector<vector<glm::vec3>> &sorPoints) {
    vector<vector<glm::vec3>> reorderedPoints;
    double pi = 2 * acos(0.0);
    reorderedPoints.resize(round(2 * pi / 0.1), vector<glm::vec3>(sorPoints[0].size()));
    int firstR = 0;
    int lastR = sorPoints[0].size() - 1;
    bool wrap = 0;
    for (int i = 0; i < sorPoints.size(); i++) {
        reorderedPoints[i].clear();
        if (i == sorPoints.size() - 1) {
            wrap = 1;
        }
        for (int j = 0; j < sorPoints[0].size(); j++) {
            if (j == firstR) {
                if (wrap) {
                    reorderedPoints[i].push_back(sorPoints[i][j]);
                    reorderedPoints[i].push_back(sorPoints[i][j + 1]);
                    reorderedPoints[i].push_back(sorPoints[0][j]);
                } else {
                    reorderedPoints[i].push_back(sorPoints[i][j]);
                    reorderedPoints[i].push_back(sorPoints[i][j + 1]);
                    reorderedPoints[i].push_back(sorPoints[i + 1][j]);
                }
            } else if (j == lastR) {
                if (wrap) {
                    reorderedPoints[i].push_back(sorPoints[i][j]);
                    reorderedPoints[i].push_back(sorPoints[0][j - 1]);
                    reorderedPoints[i].push_back(sorPoints[0][j]);
                } else {
                    reorderedPoints[i].push_back(sorPoints[i][j]);
                    reorderedPoints[i].push_back(sorPoints[i + 1][j - 1]);
                    reorderedPoints[i].push_back(sorPoints[i + 1][j]);
                }

            } else {
                if (wrap) {
                    reorderedPoints[i].push_back(sorPoints[i][j]);
                    reorderedPoints[i].push_back(sorPoints[i][j + 1]);
                    reorderedPoints[i].push_back(sorPoints[0][j]);

                    reorderedPoints[i].push_back(sorPoints[i][j]);
                    reorderedPoints[i].push_back(sorPoints[0][j - 1]);
                    reorderedPoints[i].push_back(sorPoints[0][j]);
                } else {
                    reorderedPoints[i].push_back(sorPoints[i][j]);
                    reorderedPoints[i].push_back(sorPoints[i][j + 1]);
                    reorderedPoints[i].push_back(sorPoints[i + 1][j]);

                    reorderedPoints[i].push_back(sorPoints[i][j]);
                    reorderedPoints[i].push_back(sorPoints[i + 1][j - 1]);
                    reorderedPoints[i].push_back(sorPoints[i + 1][j]);
                }
            }
        }
    }
    sorPoints = reorderedPoints;
}

void genSorCurve(CPU_Geometry curvePoints, vector<vector<glm::vec3>> &sorPoints) {
    int i = 0, j = 0;
    double pi = 2 * acos(0.0);
    sorPoints.resize(round(2 * pi / 0.1), vector<glm::vec3>(curvePoints.verts.size()));

    for (float u = 0; u < 2 * pi; u += 0.1) {
        for (int v = 0; v < curvePoints.verts.size(); v += 1) {
            sorPoints[i][j] = glm::vec3(curvePoints.verts[v].x * cosf(u), curvePoints.verts[v].y,
                                        curvePoints.verts[v].x * sinf(u));
            j++;
        }
        j = 0;
        i++;
    }
    reorderSorCurve(sorPoints);
}

glm::vec3 tsDecastlejau(vector<glm::vec3> controlPoints, int degree, float u){
    for (int i = 1; i < degree; i++) {
        for (int j = 0; j < (degree - i); j++) {
            controlPoints[j] = (1 - u) * controlPoints[j] + u * controlPoints[j + 1];
        }
    }
    if (controlPoints.size() > 0) {
        return controlPoints[0];
    }
}

//void reorderTensorSurface(vector<vector<glm::vec3>>& tensorPoints) {
//    vector<vector<glm::vec3>> reorderedPoints;
//    reorderedPoints.resize(tensorPoints.size(), vector<glm::vec3>(tensorPoints[0].size()));
//
//    int firstR = 0;
//    int lastR = tensorPoints[0].size() - 1;
//
//    for (int i = 0; i < tensorPoints.size() - 1; i++) {
//        reorderedPoints[i].clear();
//        for (int j = 0; j < tensorPoints[0].size(); j++) {
//			if (j == firstR) {
//				reorderedPoints[i].push_back(tensorPoints[i][j]);
//				reorderedPoints[i].push_back(tensorPoints[i][j + 1]);
//				reorderedPoints[i].push_back(tensorPoints[i + 1][j]);
//			
//            } else if (j == lastR) {
//                reorderedPoints[i].push_back(tensorPoints[i][j]);
//                reorderedPoints[i].push_back(tensorPoints[i + 1][j - 1]);
//                reorderedPoints[i].push_back(tensorPoints[i + 1][j]);
//            } else {
//                reorderedPoints[i].push_back(tensorPoints[i][j]);
//                reorderedPoints[i].push_back(tensorPoints[i][j + 1]);
//                reorderedPoints[i].push_back(tensorPoints[i + 1][j]);
//
//                reorderedPoints[i].push_back(tensorPoints[i][j]);
//                reorderedPoints[i].push_back(tensorPoints[i + 1][j - 1]);
//                reorderedPoints[i].push_back(tensorPoints[i + 1][j]);
//            }
//        }
//    }
//    tensorPoints = reorderedPoints;
//}
void reorderTensorSurface(vector<vector<glm::vec3>>& tensorPoints) {
	vector<vector<glm::vec3>> reorderedPoints;
	reorderedPoints.resize((tensorPoints.size() - 1) * (tensorPoints[0].size() - 1) * 2, vector<glm::vec3>(3));

	int k = 0;
	for (int i = 0; i < tensorPoints.size() - 1; i++) {
		for (int j = 0; j < tensorPoints[0].size() - 1; j++) {
			reorderedPoints[k][0] = tensorPoints[i][j];
			reorderedPoints[k][1] = tensorPoints[i][j + 1];
			reorderedPoints[k][2] = tensorPoints[i + 1][j];
			k++;

			reorderedPoints[k][0] = tensorPoints[i + 1][j];
			reorderedPoints[k][1] = tensorPoints[i][j + 1];
			reorderedPoints[k][2] = tensorPoints[i + 1][j + 1];
			k++;
		}
	}
	tensorPoints.clear();
	tensorPoints.resize(reorderedPoints.size(), vector<glm::vec3>(3));
	for (int i = 0; i < reorderedPoints.size(); i++) {
		tensorPoints[i][0] = reorderedPoints[i][0];
		tensorPoints[i][1] = reorderedPoints[i][1];
		tensorPoints[i][2] = reorderedPoints[i][2];
	}
}


////////////////////////////////////////
void genTensorSurface(vector<vector<glm::vec3>> &tensorPoints, vector<vector<glm::vec3>> &tsPoints, int d) {
    int r = 0;
    vector<vector<glm::vec3>> result;

    for (float u = 0; u < 1.f; u += 0.1f) {
        result.push_back(vector<glm::vec3>());
        int s = 0;

        for (float v = 0; v < 1.f; v += 0.1f) {
            glm::vec3 tsp = glm::vec3();
            vector<glm::vec3> R = vector<glm::vec3>();

            for (int i = 0; i < d; i++) {
                R.push_back(tsDecastlejau(tensorPoints[i], tensorPoints.size(), v));
            }
            result[r].push_back(tsDecastlejau(R, R.size(), u));
            s++;
        }
        r++;
    }
    tsPoints = result;
    reorderTensorSurface(tsPoints);
}



int getSquareIndex(CPU_Geometry square, glm::vec2 mouseCoord) {
	
    for (int i = 0; i < square.verts.size(); i++) {
        bool overlapX = (square.verts[i].x - 0.05f < mouseCoord.x) && (mouseCoord.x < square.verts[i].x + 0.05f);
        bool overlapY = (square.verts[i].y - 0.05f < mouseCoord.y) && (mouseCoord.y < square.verts[i].y + 0.05f);
        if (overlapX && overlapY) {
			
            return i;
        }
    }
	return -1;
}

std::vector<int> getTensorIndex(std::vector<std::vector<glm::vec3>> square, glm::vec2 mouseCoord) {
	std::vector<int> indexes;
	//mouseCoord = mouseCoord
	
	for (int i = 0; i < square.size(); i++) {
		for (int j = 0; j < square.size(); j++) {
			bool overlapX = (square[i][j].x - 0.05f < mouseCoord.x) && (mouseCoord.x < square[i][j].x + 0.05f);
			bool overlapY = (square[i][j].y - 0.05f < mouseCoord.y) && (mouseCoord.y < square[i][j].y + 0.05f);
			if (overlapX && overlapY) {

				//std::cout <<"mousecoords:  "<< mouseCoord.x << " , " << mouseCoord.y<<std::endl;

				indexes.push_back(i);
				indexes.push_back(j);
				return indexes;
				
			}
		}
		
	}
	indexes.push_back(-1);
	indexes.push_back(-1);
	return indexes;
	
}

/*
 * https://learnopengl.com/Getting-started/Camera - referenced
 */
void moveCamera(GLint modelLoc, GLint viewLoc, GLint projectionLoc, int view_index) {
    glm::mat4 model = glm::mat4(1.0f);
    model = glm::rotate(model, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
	glm::mat4 view;
	if (view_index == 1) {
		view = glm::lookAt(glm::vec3(5.0f, 0.0f, 0.0f), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
		view_index = 0;
	}
	else if (view_index == 2) {
		view = glm::lookAt(glm::vec3(0.0f, 0.0f, -5.0f), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
		view_index = 0;
	}
	else{
		view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
	}
    
    glm::mat4 projection = glm::perspective(glm::radians(45.0f), 800.0f / 800.0f, 0.1f, 100.0f);
	//print_mat4(view);
    glUniformMatrix4fv(modelLoc, true, 1, &model[0][0]);
    glUniformMatrix4fv(viewLoc, true, 1, &view[0][0]);
    glUniformMatrix4fv(projectionLoc, true, 1, &projection[0][0]);
}

/*
 * https://learnopengl.com/Getting-started/Camera - referenced
 */
void mouseCamera(glm::vec2 mouseCoord) {
    float xpos = mouseCoord.x;
    float ypos = mouseCoord.y;
    if (firstMouse) {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top
    lastX = xpos;
    lastY = ypos;

    float sensitivity = 0.15f; // change this value to your liking
	//float sensitivity = 1.f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    yaw += xoffset;
    pitch += yoffset;

    // make sure that when pitch is out of bounds, screen doesn't get flipped
    if (pitch > 89.0f)
        pitch = 89.0f;
    if (pitch < -89.0f)
        pitch = -89.0f;

    glm::vec3 front;
    front.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
    front.y = sin(glm::radians(pitch));
    front.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
    cameraFront = glm::normalize(front);
}


vector<vector<glm::vec3>> controlPointsGenerator( int number_of_rows) {
	vector<vector<glm::vec3>> controlPoints;
	controlPoints.resize(number_of_rows, vector<glm::vec3>(number_of_rows));
	for (int i = 0; i < number_of_rows; i++) {
		for (int j = 0; j < number_of_rows; j++) {
			controlPoints[i][j] = glm::vec3(-0.75 + i * (1.5f) / number_of_rows, -0.75 + j * (1.5f) / number_of_rows, 0.0);
		}
	}
	return controlPoints;
}

int main() {

	


	vector<vector<glm::vec3>> controlPoints;
	vector<vector<glm::vec3>> tsPoints;
	controlPoints = controlPointsGenerator(10);




    Log::debug("Starting main");

    // WINDOW
    glfwInit();
    Window window(800, 800, "CPSC 453"); // can set callbacks at construction if desired

    GLDebug::enable();
	/*if (state.scene == 7) {
		
	}*/
	
    // CALLBACKS
    shared_ptr<MyCallbacks> a3 = make_shared<MyCallbacks>(width, height);
    window.setCallbacks(a3);

	State& state = a3->getState();
	int selectedIndex = -1;
	ShaderProgram shader("shaders/test.vert", "shaders/test.frag");
	ShaderProgram heightmap_shader("shaders/heatmap.vert", "shaders/heatmap.frag");
	
	
		
	



	window.setupImGui();
    // CPU
    CPU_Geometry square;
    CPU_Geometry curve;
    CPU_Geometry sorCurve;
    CPU_Geometry tensorControl;
    CPU_Geometry tensorSurface;
	CPU_Geometry_height heightMap;

    // GPU
    GPU_Geometry pointsGPUGeom;
    GPU_Geometry linesGPUGeom;
    GPU_Geometry curveGPUGeom;
    GPU_Geometry sorGPUGeom;
    GPU_Geometry tensorCGUPGeom;
    GPU_Geometry tensorGPUGeom;
	GPU_Geometry_height heightMapGPU;


	std::vector<float> heightValues;
	float height_norm = 0.0f;
	
    glPointSize(10.0f);
    // RENDER LOOP
    while (!window.shouldClose()) {
		float minHeight = 0.f; float maxHeight = 0.f;
		heightMap.verts.clear();
		heightMap.cols.clear();
		heightMap.heights.clear();
		for (glm::vec3 point : tensorSurface.verts) {
		
			heightMap.heights.push_back(point.z);
			//std::cout << point.z << std::endl;
		}
		/*int a;
		std::cin >> a;*/
		if (heightMap.heights.size() == 0) {
			minHeight = 0.f;
			maxHeight = 0.f;
		}
		else {
			auto minmax = std::minmax_element(heightMap.heights.begin(), heightMap.heights.end());
			minHeight = *minmax.first;
			maxHeight = *minmax.second;

		}
		
		for (float point : heightMap.heights) {
			glm::vec3 color;
			if (abs(maxHeight - minHeight) < 0.001) {
				height_norm = 0.0;
			}
			else {
				height_norm = (point - minHeight) / (maxHeight - minHeight);
				/*if (height_norm > 0.0f) {
					std::cout << height_norm << std::endl;
				}*/
				//std::cout << height_norm << std::endl;
			}
			if (point < threshold2) {
				color = mix(color1, color2, (point - threshold1) / (threshold2 - threshold1));
			}
			else if (point < threshold3) {
				color = mix(color3, color2, (point - threshold2) / (threshold3 - threshold2));
			}
			else if (point < threshold4) {
				color = mix(color4, color3, (point - threshold3) / (threshold4 - threshold3));
			}
			else {
				color = mix(color4, color5, (point - threshold4) / (threshold5 - threshold4));
			}
			//color = glm::vec3(1 - height_norm, 0.f, height_norm);
			heightMap.cols.push_back(color);
		}
		
		

        glfwPollEvents();

		// Start the Dear ImGui frame
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		ImGui::Begin("Handles");
		ImGui::Text("blah blah blah");
		if (ImGui::Button("Normal view")) {
			view_index = 0;
		}
		if (ImGui::Button("View 1")) {
			view_index = 1;
		}
		if (ImGui::Button("View 2")) {
			view_index = 2;
		}

		// Add input fields for the light position
		ImGui::InputFloat3("Position", &lightPosition[0]);


		ImGui::End();

		ImGui::Render();




        GLint modelLoc = glGetUniformLocation(shader.getProgram(), "model");
        GLint viewLoc = glGetUniformLocation(shader.getProgram(), "view");
        GLint projectionLoc = glGetUniformLocation(shader.getProgram(), "projection");
		GLint lightpositionLoc = glGetUniformLocation(shader.getProgram(), "lightposition");


		GLint modelLoc1 = glGetUniformLocation(heightmap_shader.getProgram(), "model1");
		GLint viewLoc1 = glGetUniformLocation(heightmap_shader.getProgram(), "view1");
		GLint projectionLoc1 = glGetUniformLocation(heightmap_shader.getProgram(), "projection1");
		//GLint heightValuesLoc = glGetUniformLocation(heightmap_shader.getProgram(), "heightValues");
		// Get the location of the uniform variables in the shaders
		//GLint maxHeightLoc = glGetUniformLocation(heightmap_shader.getProgram(), "maxHeight");
		//GLint minHeightLoc = glGetUniformLocation(heightmap_shader.getProgram(), "minHeight");
        //GLint modeLoc = glGetUniformLocation(shader.getProgram(), "mode");

		glUniform3f(lightpositionLoc, lightPosition[0], lightPosition[1], lightPosition[2]);
		//glUniform1fv(heightValuesLoc, heightValues.size(), &heightValues[0]);
		//glUniform1f(maxHeightLoc, maxHeight);
		//glUniform1f(minHeightLoc, minHeight);

	        if (a3->checkStateChanged()) {
            if (state.leftMousePressed && selectedIndex == -1 && !state.mode) {
                selectedIndex = getSquareIndex(square, state.mouseCoord);
                if (selectedIndex == -1) {
                    square.verts.push_back(glm::vec3(state.mouseDownCoord, 0.f));
                    square.cols.resize(square.verts.size(), glm::vec3{1.0, 0.0, 0.0});
                }
            } else if (state.leftMouseHeld && !state.mode) {
                if (selectedIndex == -1) {
                    square.verts.back() = glm::vec3(state.mouseCoord, 0.f);
                    square.cols.resize(square.verts.size(), glm::vec3{1.0, 0.0, 0.0});
                } else {
                    square.verts[selectedIndex] = glm::vec3(state.mouseCoord, 0.f);
                }
            } else if (state.leftMouseReleased && !state.mode) {
                selectedIndex = -1;
            } else if (state.rightMousePressed && selectedIndex == -1 && !state.mode) {
                selectedIndex = getSquareIndex(square, state.mouseCoord);
                if (selectedIndex != -1) {
                    square.verts.erase(square.verts.begin() + selectedIndex);
                    selectedIndex = -1;
                }
            }




			else if (state.leftMouseHeld && state.mode && state.pups && state.wire && !state.controlDown && !state.controlUp) {
				/*if (selectedIndex == -1) {
					tensorControl.verts.back() = glm::vec3(state.mouseCoord, 0.f);
					tensorControl.cols.resize(tensorControl.verts.size(), glm::vec3{ 1.0, 0.0, 0.0 });
				}*/
				std::vector<int> index = getTensorIndex(controlPoints, state.mouseCoord);
				if(index[0] != -1) {
					controlPoints[index[0]][index[1]] = glm::vec3(state.mouseCoord, controlPoints[index[0]][index[1]].z);
					//std::cout <<"control Points: "<< controlPoints[index[0]][index[1]].x << " , " << controlPoints[index[0]][index[1]].y << " , " << controlPoints[index[0]][index[1]].z << std::endl;
					
				}
			}
			else if (state.leftMouseHeld && state.mode && state.pups && state.wire && !state.controlDown && state.controlUp) {
				/*if (selectedIndex == -1) {
					tensorControl.verts.back() = glm::vec3(state.mouseCoord, 0.f);
					tensorControl.cols.resize(tensorControl.verts.size(), glm::vec3{ 1.0, 0.0, 0.0 });
				}*/
				std::vector<int> index = getTensorIndex(controlPoints, state.mouseCoord);
				if (index[0] != -1) {
					controlPoints[index[0]][index[1]] = glm::vec3(state.mouseCoord, controlPoints[index[0]][index[1]].z + 0.01);
					//std::cout << "control Points: " << controlPoints[index[0]][index[1]].x << " , " << controlPoints[index[0]][index[1]].y << " , " << controlPoints[index[0]][index[1]].z << std::endl;

				}
				state.controlUp = false;
			}
			else if (state.leftMouseHeld && state.mode && state.pups && state.wire && state.controlDown && !state.controlUp) {
				/*if (selectedIndex == -1) {
					tensorControl.verts.back() = glm::vec3(state.mouseCoord, 0.f);
					tensorControl.cols.resize(tensorControl.verts.size(), glm::vec3{ 1.0, 0.0, 0.0 });
				}*/
				std::vector<int> index = getTensorIndex(controlPoints, state.mouseCoord);
				if (index[0] != -1) {
					controlPoints[index[0]][index[1]] = glm::vec3(state.mouseCoord, controlPoints[index[0]][index[1]].z - 0.01);
					//std::cout << "control Points: " << controlPoints[index[0]][index[1]].x << " , " << controlPoints[index[0]][index[1]].y << " , " << controlPoints[index[0]][index[1]].z << std::endl;

				}
				state.controlDown = false;
			}
			else if (state.leftMouseReleased && state.mode && state.pups && state.wire) {
				selectedIndex = -1;
			}




			else if (state.scene == -1) {
                square.verts.clear();
                square.cols.clear();
                curve.verts.clear();
                curve.cols.clear();
                sorCurve.verts.clear();
                sorCurve.cols.clear();
				tensorControl.verts.clear();
				tensorControl.cols.clear();
				tensorControl.normals.clear();
				tensorSurface.verts.clear();
				tensorSurface.cols.clear();
				tensorSurface.normals.clear();
                a3->resetScene();
            }
    //        // pups
    //        if (state.scene == 1) {
    //            vector<glm::vec3> controlPoints;
    //            vector<glm::vec3> PupsPoints;

    //            controlPoints.assign(square.verts.begin(), square.verts.end());
				//for (float u = 0; u < 1; u += 0.01) {
				//	
				//	PupsPoints.push_back(pups(controlPoints, u));

				//}
				///*
				//for (int j = 0; j < PupsPoints.size(); j += 1) {

				//	std::cout << "point : (" << PupsPoints[j].x << ", " << PupsPoints[j].y << ", " << PupsPoints[j].z << ")\n";
				//}
				//*/
				////std::cout << PupsPoints.size() << "\n";
    //            sorCurve.verts.clear();
    //            sorCurve.cols.clear();
    //            curve.verts.assign(PupsPoints.begin(), PupsPoints.end());
    //            curve.cols.clear();
    //            curve.cols.resize(curve.verts.size(), glm::vec3{2.0f, 0.5f, 1.0f});
    //        }
            
            
            
			if (state.scene == 7 || state.scene ==2) {
				tensorSurface.verts.clear();
				

				vector<glm::vec3> tempSpace;
				
				float ratio = 20;
				for (float u = 0; u <= 1; u += (1/(ratio-1))) {
					for (float v = 0; v <= 1; v += (1 / (ratio-1))) {
						tempSpace.push_back(pups_surface(controlPoints, u, v));
					}
				}


				tensorControl.verts.clear();
				tensorControl.cols.clear();
				tensorControl.normals.clear();

				//tensorControl.verts = tempSpace;

				for (int i = 0; i < controlPoints.size(); i++) {
					for (int j = 0; j < controlPoints[0].size(); j++) {
						tensorControl.verts.push_back(controlPoints[i][j]);
					}
				}

				/////////////////////////////////////////////

				
				//vector<glm::vec3> surface_grid_pts;
				vector<vector<glm::vec3>> surface_grid;
				surface_grid.resize(int(ratio), vector<glm::vec3>(int(ratio)));

				for (int i = 0; i < tempSpace.size(); i++)
				{
					surface_grid[i / ratio][i % int(ratio)] = tempSpace[i];
				}

				reorderTensorSurface(surface_grid);
				//reorder_points(surface_grid);
				for (int i = 0; i < surface_grid.size(); i++) {
					for (int j = 0; j < surface_grid[i].size(); j++) {
						tensorSurface.verts.push_back(surface_grid[i][j]);
					}
				}

				
				tensorControl.cols.resize(tensorControl.verts.size(), glm::vec3{1.0f, 0.0f, 0.0f}); // tempSpace.size() tensorControl.verts.size() 
				tensorControl.normals.resize(tensorControl.verts.size(), glm::vec3{ 0.0f, 0.0f, 1.0f });
				/*for (int i = tensorControl.verts.size() - 1; i > tensorControl.verts.size() - 17; i--)
				{
					tensorControl.cols[i] = glm::vec3{ 1.0f, 0.0f, 1.0f };
				}*/
				
				tensorSurface.cols.clear();
				tensorSurface.normals.clear();
				tensorSurface.cols.resize(tensorSurface.verts.size(), glm::vec3{ 2.0f, 0.5f, 1.0f });
				tensorSurface.normals = get_normals(tensorSurface.verts);
			}
            // 3d
            if (state.mode) {
				if (state.cameraReset && !state.wire) {
					cameraPos = glm::vec3(0.0f, 0.0f, 0.0f);
					cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
					cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);

				}
                if (state.cameraForward && !state.wire) {
                    cameraPos += cameraSpeed * cameraFront;
                }
                if (state.cameraBackward && !state.wire) {
                    cameraPos -= cameraSpeed * cameraFront;
                }
                if (state.cameraLeft && !state.wire) {
                    cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
                }
                if (state.cameraRight && !state.wire) {
                    cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
                }
                if (state.cameraUp && !state.wire) {
                    cameraPos += cameraSpeed * cameraUp;
                }
                if (state.cameraDown && !state.wire) {
                    cameraPos -= cameraSpeed * cameraUp;
                }
                if (state.leftMouseHeld && !state.wire) {
                    mouseCamera(glm::vec2(a3->getMousePos()));
                }
                //if (state.revolve && (state.scene == 1 || state.scene == 2)) {
                //    vector<vector<glm::vec3>> sorPoints;
                //    genSorCurve(curve, sorPoints);
                //    for (int i = 0; i < sorPoints.size(); i++) {
                //        for (int j = 0; j < sorPoints[0].size(); j++) {
                //            sorCurve.verts.push_back(sorPoints[i][j]);
                //        }
                //    }
                //    /*if (state.scene == 1) {
                //        sorCurve.cols.clear();
                //        sorCurve.cols.resize(sorCurve.verts.size(), glm::vec3{2.0f, 0.5f, 1.0f});
                //    } else if (state.scene == 2) {
                //        sorCurve.cols.clear();
                //        sorCurve.cols.resize(sorCurve.verts.size(), glm::vec3{1.0f, 0.5f, 0.0f});
                //    }*/
                //}
                if (!state.revolve || !state.mode) {
                    sorCurve.verts.clear();
                    sorCurve.cols.clear();
                }
                /**
                 * https://stackoverflow.com/questions/137629/how-do-you-render-primitives-as-wireframes-in-opengl - referenced
                 */
                if (state.wire) {
                    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                }
                if (!state.wire) {
                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                }
            }
            a3->stateHandled();
        }
        //glUniform1i(modeLoc, state.mode);

		heightMap.verts = tensorSurface.verts;
		
        
		if (state.scene == 2) {
			heightmap_shader.use();
		}
		else {
			shader.use();
		}

		if (state.scene == 2) {
			moveCamera(modelLoc1, viewLoc1, projectionLoc1, view_index);
		}
		else {
			moveCamera(modelLoc, viewLoc, projectionLoc, view_index);
		}
        updateGPUGeometry(pointsGPUGeom, square);
        updateGPUGeometry(linesGPUGeom, square);
        updateGPUGeometry(curveGPUGeom, curve);
        updateGPUGeometry(sorGPUGeom, sorCurve);
        updateGPUGeometry(tensorCGUPGeom, tensorControl);
        updateGPUGeometry(tensorGPUGeom, tensorSurface);
		updateHeightGPUGeometry(heightMapGPU, heightMap);
        glEnable(GL_LINE_SMOOTH);
        glEnable(GL_FRAMEBUFFER_SRGB);
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


		

        if (!state.mode) {
            linesGPUGeom.bind();
            glDrawArrays(GL_LINE_STRIP, 0, GLsizei(square.verts.size()));

            pointsGPUGeom.bind();
            glDrawArrays(GL_POINTS, 0, GLsizei(square.verts.size()));
        }

        if (state.scene == 1 ){
            curveGPUGeom.bind();
            glDrawArrays(GL_LINE_STRIP, 0, GLsizei(curve.verts.size()));

            sorGPUGeom.bind();
            glDrawArrays(GL_TRIANGLES, 0, GLsizei(sorCurve.verts.size()));
        } else if (state.mode){
			if (state.scene == 7) {
				tensorCGUPGeom.bind();
				glDrawArrays(GL_POINTS, 0, GLsizei(tensorControl.verts.size()));

				tensorGPUGeom.bind();
				glDrawArrays(GL_TRIANGLES, 0, GLsizei(tensorSurface.verts.size()));
			}
			if (state.scene == 2) {
				tensorCGUPGeom.bind();
				glDrawArrays(GL_POINTS, 0, GLsizei(tensorControl.verts.size()));

				heightMapGPU.bind();
				glDrawArrays(GL_TRIANGLES, 0, GLsizei(tensorSurface.verts.size()));
			}
        }

        glDisable(GL_FRAMEBUFFER_SRGB); // disable sRGB for things like imgui

		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        window.swapBuffers();
    }

	// Cleanup
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();


    glfwTerminate();
    return 0;
}
