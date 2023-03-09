#include "Geometry.h"

#include <utility>


GPU_Geometry::GPU_Geometry()
	: vao()
	, vertBuffer(0, 3, GL_FLOAT)
	, colorsBuffer(1, 3, GL_FLOAT)
	, normalsBuffer(2, 3, GL_FLOAT)
{}


void GPU_Geometry::setVerts(const std::vector<glm::vec3>& verts) {
	vertBuffer.uploadData(sizeof(glm::vec3) * verts.size(), verts.data(), GL_STATIC_DRAW);
}


void GPU_Geometry::setCols(const std::vector<glm::vec3>& cols) {
	colorsBuffer.uploadData(sizeof(glm::vec3) * cols.size(), cols.data(), GL_STATIC_DRAW);
}

void GPU_Geometry::setNormals(const std::vector<glm::vec3>& norms) {
	normalsBuffer.uploadData(sizeof(glm::vec3) * norms.size(), norms.data(), GL_STATIC_DRAW);
}

GPU_Geometry_height::GPU_Geometry_height()
	: vao()
	, vertBuffer(0, 3, GL_FLOAT)
	, heightBuffer(1, 1, GL_FLOAT)
	, colorsBuffer(2, 3, GL_FLOAT)
{}


void GPU_Geometry_height::setVerts(const std::vector<glm::vec3>& verts) {
	vertBuffer.uploadData(sizeof(glm::vec3) * verts.size(), verts.data(), GL_STATIC_DRAW);
}

void GPU_Geometry_height::setCols(const std::vector<glm::vec3>& cols) {
	colorsBuffer.uploadData(sizeof(glm::vec3) * cols.size(), cols.data(), GL_STATIC_DRAW);
}
void GPU_Geometry_height::setHeight(const std::vector<float>& heights) {
	heightBuffer.uploadData(sizeof(float) * heights.size(), heights.data(), GL_STATIC_DRAW);
}

