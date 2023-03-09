#version 330 core
    in float aheight;
	//in float OmaxHeight;
	//in float OminHeight;
	in vec3 fragColor;


    out vec4 FragColor;
    //uniform float heightValues[512 * 512]; // Assume 'heightValues' is an array of height values for each point
	



	// Define the color thresholds and corresponding colors
	//const vec3 color1 = vec3(0.0, 0.0, 0.5);   // dark blue
	//const vec3 color2 = vec3(0.0, 0.5, 1.0);   // light blue
	//const vec3 color3 = vec3(0.9, 0.9, 0.0);   // yellow
	//const vec3 color4 = vec3(0.5, 0.0, 0.0);   // dark red
	//const vec3 color5 = vec3(1.0, 0.0, 0.0);   // bright red



    void main()
    {
		
		// Normalize the height value to the range [0, 1]
		//float height_norm;
		//if (OmaxHeight == OminHeight) {
		//	 height_norm = 0.0;
		//} else {
		//	 height_norm = (aheight - OminHeight) / (OmaxHeight - OminHeight);
		//}
		//// Define the height thresholds for each color
		//float threshold1 = 0.0;
		//float threshold2 = 0.1;
		//float threshold3 = 0.2;
		//float threshold4 = 0.3;
		//float threshold5 = 1.0;
		// Interpolate between colors based on the height value
		//vec3 color;
		//if (height_norm < threshold2) {
		//    color = color1;
		//} else if (height_norm < threshold3) {
		//    color = color2;
		//} else if (height_norm < threshold4) {
		//    color = color3;
		//} else {
		//    color = color4;
		//}
		//color = mix(color4, color5, (height_norm-threshold4)/(threshold5-threshold4));
		//FragColor = vec4(color, 1.0);
		//FragColor = vec4(1- height_norm, 0.f, height_norm, 1.0);
		//FragColor = vec4(vec3(0.5, 0.0,  1.0), 1.0);
		// Get the height value for the current point
		//float height = heightValues[int(gl_FragCoord.y * 512 + gl_FragCoord.x)];

		// Normalize the height value to the range [0, 1]
		//float height_norm = (height - minHeight) / (maxHeight - minHeight);

		// Calculate the color based on the height value
		//vec3 color = vec3(height_norm, 0.0, 1-height_norm);

		//FragColor = vec4(color, 1.0);
		//FragColor = vec4(1.f,0.f,0.f, 1.0);



        //float colorIndex = height * 4.0;
        //vec3 color = texture(colormap, colorIndex).rgb;
		//color = height
        //FragColor = vec4(color, 1.0f);



		FragColor = vec4(fragColor, 1.0);
    }
