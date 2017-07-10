function [ serial_obj ] = UVConn(port)
    %This code creates the serial object. The 'fopen' line is connecting to
    %the UV lamp. The code combines serial object creation and connection.
    %Alternative serial connection code would look like: 
    % fprintf(serial_obj,['CONN18',char(13)]);

    serial_obj = serial(port,'BAUD',19200,'DataBits',8,'StopBits',1,'Parity','None','Terminator','');
    fopen(serial_obj); 
    
end

