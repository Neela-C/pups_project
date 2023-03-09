#version 330 core
    layout (location = 0) in vec3 aPos;
    layout (location = 1) in float aHeight;
	layout (location = 2) in vec3 color;

    uniform mat4 model1;
    uniform mat4 view1;
    uniform mat4 projection1;
	//uniform float maxHeight; // The maximum height value in the array
	//uniform float minHeight; // The minimum height value in the array
    out float aheight;
	//out float OmaxHeight;
	//out float OminHeight;


	out vec3 fragColor;


    void main()
    {
		mat4 M =model1 * view1 * projection1 ;
        gl_Position =  vec4(aPos, 1.0) * M;
	    //mat4 M = projection1 * view1 * model1 ;
	    //gl_Position = M* vec4(aPos, 1.0);
        aheight = aHeight;
		//OmaxHeight = maxHeight;
		//OminHeight = minHeight;
		fragColor = color;
    }
