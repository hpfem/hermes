var popupHandle;

function popup(myPopup)
{
  var win = "videos.html" + myPopup;
  // always close the old one , so only one at a time is open
  if(popupHandle || popupHandle!=null)
  {
    if (!popupHandle.closed) popupHandle.close();
  }
  popupHandle=null;
   if(true)
 {
  // open the popup
  popupHandle = window.open( win, "myWindow", 
"height = 550, width = 900, resizable = 0" )
  return popupHandle;
 }
}

function winclose()
{
  if (window.popupHandle!=null && !window.popupHandle.closed)
  {
    window.popupHandle.close();
    }
}

function doNothing(){} // does nothing but required by JavaScript in this case

