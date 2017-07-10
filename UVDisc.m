function UVDisc(serial_obj)
    % Similar to UVConn, this code uses 'fclose' to disconnect from the
    % serial object. Alternative serial disconnection code would look like: 
    % fprintf(serial_obj,['DCONE1',char(13)]);

    fclose(serial_obj);
    
    delete(serial_obj);  %deletes the serial object
    clear serial_obj;
end