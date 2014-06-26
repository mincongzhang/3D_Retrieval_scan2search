【OpenGL】三维场景下用鼠标在window上绘图
Software: Microsoft visual studio 2012 MFC

Library: opengl

首先，自定义变量

//Transform from screen coordination to OpenGL coordination
GLint viewport[4]; 
GLdouble modelview[16]; 
GLdouble projection[16]; 
GLfloat winX, winY, winZ; 
GLdouble posX, posY, posZ;

然后获得转换矩阵

glGetIntegerv(GL_VIEWPORT, viewport);
glGetDoublev(GL_MODELVIEW_MATRIX, modelview); 
glGetDoublev(GL_PROJECTION_MATRIX, projection);

最后得到对应的OpenGL坐标

 winX = (float)x; 
 winY = viewport[3] - (float)y;
 glReadPixels((int)winX, (int)winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ); 
 gluUnProject(winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ); 

注：(x, y)是屏幕坐标，(winX, winY, winZ)是视景体坐标及深度坐标，(posX, posY, posZ)是OpenGL坐标。

此方法需在glViewport(0, 0, screenWidth, screenHeight)情况下，screenWidth、screenHeight分别是客户区的宽和高，视口左下角坐标恰好是（0，0），并且未经过任何模型变换。

方法已验证。

若非此情况，参照：http://blog.csdn.net/abcdef8c/article/details/6716737